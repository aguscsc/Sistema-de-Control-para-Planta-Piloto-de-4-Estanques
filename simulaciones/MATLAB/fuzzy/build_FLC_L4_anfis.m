%% build_FLC_L4_anfis.m
% Level loop LIC4 (upper tank 4) - cascade on FIC1
%
% Tank 4 (in % units: H = level % of h_max, F = pump 1 flow % of q_max):
%   dH/dt = beta*F - alpha*sqrt(H)
%
% Controller (positional FF + accumulated FB):
%   F_ff  = cFF*sqrt(r)                      static feedforward (exact steady state)
%   I_k   = I_{k-1} + dFB_k                  FB accumulator (conditional integration)
%   F_sp  = sat(F_ff + I_k, 0, 100)          flow setpoint for FIC1
%
% FB FLC (same method as the flow loops):
%   s   = e + Ti*de                          predicted error if flow were frozen
%   dFB = Kdu*sat(s/delta)                   -> Dahlin PI near setpoint, rate-capped far away
%   e   = r - H,  de = -(H_k - H_{k-1})/Ts   (P on measurement)
%
% Why positional FF: a velocity-form FF increment is lost whenever the
% total setpoint is clamped (0/100 %), leaving an offset the FB must undo.
% Keeping FF absolute and only accumulating the FB avoids that.

clear; clc; close all;

%% ================= USER SETTINGS =================
nMF    = 3;        % triangles per input
nOut   = 0;        % shared singletons after training (0 = keep one per rule)
Hd     = 30;       % [%] design operating point for Ti
lambdaL = 20;      % [s] desired closed-loop time constant  <-- design knob

%% ================= Plant =================
g = 9.81; A4 = 0.159; a4 = 4.1820e-04; gamma1 = 0.455;
q_max = 0.0017; h_max = 0.6; k1 = q_max;

T.beta  = (1-gamma1)*k1/(A4*h_max);                    % (%/s) per % flow
T.alpha = 100*a4*sqrt(2*g*h_max/100)/(A4*h_max);       % %/s per sqrt(%)
T.cFF   = T.alpha/T.beta;                              % F_ss = cFF*sqrt(H)
T.Hmax  = (100/T.cFF)^2;                               % max level at F = 100 %

fprintf('===== LIC4 =====\n');
fprintf('beta = %.5f   alpha = %.4f   cFF = %.3f\n', T.beta, T.alpha, T.cFF);
fprintf('Max steady level with pump at 100%% : %.1f %%\n', T.Hmax);

% Linearized at the design point
T.tau_d = 2*sqrt(Hd)/T.alpha;
T.K_d   = 2*sqrt(Hd)/T.cFF;
fprintf('At H = %g %% : tau = %.1f s   K = %.3f   (K/tau = beta)\n', Hd, T.tau_d, T.K_d);

%% ================= Controller parameters =================
P.Ts     = 0.1;
P.rate   = 1.0;                         % [%/s] flow can't move faster than this
P.lambda = lambdaL;

% Dahlin PI on the linearized tank at Hd
P.a  = exp(-P.Ts/T.tau_d);
P.q  = exp(-P.Ts/P.lambda);
P.Kc = (1 - P.q)/(T.K_d*(1 - P.a));     % % flow per % level (~ 1/(beta*lambda))
P.Ti = P.Ts/(1 - P.a);

% Target law
P.Kdu   = P.rate*P.Ts;                  % max flow-setpoint move per sample
P.delta = P.Kdu*P.Ti/(P.Kc*P.Ts);       % [% level]

% Universes
P.DEmax = 0.3;                          % [%/s] realistic level rate (beta*30 % flow excess)
P.Emax  = P.Ti*P.DEmax + 2*P.delta;     % [% level]
P.Emax  = min(P.Emax, 100);
P.Se    = 1/P.Emax;
P.Sde   = 1/P.DEmax;

fprintf('--- Dahlin PI ---\n');
fprintf('Kc = %.3f %%flow/%%level   Ti = %.1f s\n', P.Kc, P.Ti);
fprintf('--- Target law ---\n');
fprintf('Kdu = %.2f %%/sample   delta = %.2f %%\n', P.Kdu, P.delta);
fprintf('Emax = %.1f %% -> Se = %.4f   DEmax = %.2f %%/s -> Sde = %.3f (Simulink gain on -dH: %.2f)\n', ...
        P.Emax, P.Se, P.DEmax, P.Sde, P.Sde/P.Ts);

%% ================= Training data =================
target = @(en, den) max(-1, min(1, (en*P.Emax + P.Ti*den*P.DEmax)/P.delta));
gg = linspace(-1, 1, 61);
[EN, DEN] = meshgrid(gg);
Xtr = [EN(:) DEN(:)];  Ytr = target(Xtr(:,1), Xtr(:,2));
rng(1);
Xva = 2*rand(1500,2) - 1;  Yva = target(Xva(:,1), Xva(:,2));

%% ================= ANFIS =================
gopt = genfisOptions('GridPartition');
gopt.NumMembershipFunctions       = [nMF nMF];
gopt.InputMembershipFunctionType  = 'trimf';
gopt.OutputMembershipFunctionType = 'constant';
fis0 = genfis(Xtr, Ytr, gopt);

aopt = anfisOptions('InitialFIS', fis0, 'EpochNumber', 200, ...
    'ValidationData', [Xva Yva], 'DisplayANFISInformation', 0, ...
    'DisplayErrorValues', 0, 'DisplayStepSize', 0, 'DisplayFinalResults', 1);
[fisT, trErr, ~, fisV, vaErr] = anfis([Xtr Ytr], aopt);

opt  = evalfisOptions('OutOfRangeInputValueMessage','none','NoRuleFiredMessage','none');
rmse = @(f) sqrt(mean((evalfis(f, Xtr, opt) - Ytr).^2));
cand = {fis0, fisT, fisV};  cname = ["initial (genfis)" "final epoch" "min validation"];
r = cellfun(rmse, cand);
fprintf('--- evalfis RMSE of each candidate ---\n');
for i = 1:3, fprintf('%-18s %.4f\n', cname(i), r(i)); end
[~, ib] = min(r);  fis = cand{ib};  fprintf('Using: %s\n', cname(ib));

name = 'FLC_L4';
fis.Name = name;  fis.Inputs(1).Name = 'e';  fis.Inputs(2).Name = 'de';  fis.Outputs(1).Name = 'du';

nR = numel(fis.Rules);
cTab = nan(nMF);
for k = 1:nR
    ant = fis.Rules(k).Antecedent;
    cTab(ant(1),ant(2)) = fis.Outputs(1).MembershipFunctions(fis.Rules(k).Consequent).Parameters;
end
fprintf('--- Trained singletons (rows e = N..P, cols de = N..P) ---\n');  disp(round(cTab,3));

%% ================= Optional: reduce to nOut shared singletons =================
if nOut > 0
    c = zeros(nR,1);
    for k = 1:nR, c(k) = fis.Outputs(1).MembershipFunctions(fis.Rules(k).Consequent).Parameters; end
    lev0 = linspace(min(c), max(c), nOut);
    [~, L] = min(abs(c - lev0), [], 2);

    W = zeros(size(Xtr,1), nR);
    for k = 1:size(Xtr,1)
        [~,~,~,~,rf] = evalfis(fis, Xtr(k,:), opt);
        W(k,:) = rf(:)' / max(sum(rf), eps);
    end
    A = zeros(size(Xtr,1), nOut);
    for j = 1:nOut, A(:,j) = sum(W(:, L == j), 2); end
    used = any(A,1);  v = lev0;  v(used) = A(:,used) \ Ytr;

    if nOut == 5, outNames = ["NB" "NS" "ZE" "PS" "PB"];
    elseif nOut == 3, outNames = ["N" "Z" "P"];
    else, outNames = "o" + string(1:nOut); end

    fisR = sugfis('Name', name, 'AndMethod','prod', 'DefuzzificationMethod','wtaver');
    for i = 1:2
        fisR = addInput(fisR, [-1 1], 'Name', fis.Inputs(i).Name);
        for m = 1:numel(fis.Inputs(i).MembershipFunctions)
            mf = fis.Inputs(i).MembershipFunctions(m);
            fisR = addMF(fisR, fis.Inputs(i).Name, mf.Type, mf.Parameters, 'Name', mf.Name);
        end
    end
    fisR = addOutput(fisR, [min(-1,min(v)) max(1,max(v))], 'Name', 'du');
    for j = 1:nOut, fisR = addMF(fisR, 'du', 'constant', v(j), 'Name', outNames(j)); end
    ruleList = zeros(nR,5);
    for k = 1:nR, ruleList(k,:) = [fis.Rules(k).Antecedent, L(k), 1, 1]; end
    fisR = addRule(fisR, ruleList);

    fprintf('--- Reduction %d -> %d shared singletons ---\n', nR, nOut);
    fprintf('LS fit RMSE: %.4f\n', sqrt(mean((A*v(:) - Ytr).^2)));
    fprintf('Singletons : %s\n', mat2str(v,3));
    fis = fisR;
end

%% ================= Checks =================
DU = reshape(evalfis(fis, Xtr, opt), size(EN));  TG = reshape(Ytr, size(EN));
fprintf('--- Fit ---\n');
fprintf('RMSE train = %.4f   max |err| = %.4f\n', sqrt(mean((DU(:)-TG(:)).^2)), max(abs(DU(:)-TG(:))));
h = 1e-3;
ge = (evalfis(fis,[h 0],opt) - evalfis(fis,[-h 0],opt))/(2*h);
gd = (evalfis(fis,[0 h],opt) - evalfis(fis,[0 -h],opt))/(2*h);
fprintf('Origin slope e : %.2f (target %.2f)   de : %.2f (target %.2f)\n', ...
        ge, P.Emax/P.delta, gd, P.Ti*P.DEmax/P.delta);

figure('Name',[name ' surface'],'Position',[100 100 1200 380]);
subplot(1,3,1); surf(EN,DEN,TG,'EdgeColor','none'); view(-35,30); title('Target law'); xlabel('e_n'); ylabel('de_n');
subplot(1,3,2); surf(EN,DEN,DU,'EdgeColor','none'); view(-35,30); title('FIS'); xlabel('e_n'); ylabel('de_n');
subplot(1,3,3); contourf(EN,DEN,DU-TG,20,'LineColor','none'); colorbar; title('Error'); axis square;

%% ================= Closed-loop test (nonlinear tank + approximate flow loop) =================
% Flow loop approximated as: command moves toward F_sp at <= 1 %/s, then pump lag
tauP = 1.3776;  aF = exp(-P.Ts/tauP);
Tsim = 900;  N = round(Tsim/P.Ts);  t = (0:N-1)*P.Ts;
r  = 20*ones(1,N);  r(t >= 50) = 30;              % setpoint step 20 -> 30 %
dF = zeros(1,N);    dF(t >= 500) = -5;            % disturbance: -5 % flow reaching tank (model error)

cases = ["FLC + FF", "PI + FF", "FF only"];
Hs = zeros(3,N);  Fs = zeros(3,N);
for cs = 1:3
    H = 20*ones(1,N);  F0 = T.cFF*sqrt(20);
    F = F0;  Fc = F0;  I = 0;  Hprev = 20;
    for k = 2:N
        % tank
        H(k) = H(k-1) + P.Ts*(T.beta*max(F + dF(k-1), 0) - T.alpha*sqrt(max(H(k-1),0)));
        H(k) = min(max(H(k), 0), 100);

        % outer controller
        Fff = T.cFF*sqrt(max(r(k),0));
        e   = r(k) - H(k);
        de  = -(H(k) - Hprev)/P.Ts;   Hprev = H(k);
        switch cs
            case 1
                en  = max(-1, min(1, P.Se*e));
                den = max(-1, min(1, P.Sde*de));
                dFB = P.Kdu*max(-1, min(1, evalfis(fis,[en den],opt)));
            case 2
                dFB = -P.Kc*(-de*P.Ts) + P.Kc*(1-P.a)*e;     % P on measurement
                dFB = max(-P.Kdu, min(P.Kdu, dFB));
            case 3
                dFB = 0;
        end
        Fsp_try = Fff + I + dFB;
        if ~((Fsp_try > 100 && dFB > 0) || (Fsp_try < 0 && dFB < 0))
            I = I + dFB;                                     % conditional integration
        end
        Fsp = min(max(Fff + I, 0), 100);

        % inner flow loop (approx.)
        Fc = Fc + max(-P.rate*P.Ts, min(P.rate*P.Ts, Fsp - Fc));
        F  = aF*F + (1-aF)*Fc;
        Fs(cs,k) = F;
    end
    Hs(cs,:) = H;  Fs(cs,1) = F0;
end

figure('Name','LIC4 closed loop');
subplot(2,1,1); plot(t, r, 'k--'); hold on; plot(t, Hs, 'LineWidth', 1.2);
legend(['setpoint' cases], 'Location','southeast'); ylabel('H_4 [%]'); grid on;
title('Step 20 \rightarrow 30 % at t = 50 s, flow disturbance -5 % at t = 500 s');
subplot(2,1,2); plot(t, Fs, 'LineWidth', 1.2); legend(cases); ylabel('F_1 [%]'); xlabel('t [s]'); grid on;

fprintf('--- Step 20->30 (2 %% band = 0.2 %% level) ---\n');
for cs = 1:3
    seg = t < 500;
    fprintf('%-9s overshoot %.2f %%   settling %.0f s   offset after dist. %.2f %%\n', cases(cs), ...
        max(Hs(cs,seg)) - 30, settleTime(Hs(cs,seg), t(seg), 30, 0.2, 50), 30 - Hs(cs,end));
end

%% ================= Save =================
writeFIS(fis, name);
FLC_L4 = fis;  L4 = P;  TK4 = T;
save('FLC_L4_params.mat', 'FLC_L4', 'L4', 'TK4');

% Simulink (LIC4, cascade on FIC1):
%  Tank 4:   F1 -> [beta] -> (+) -> Integrator (IC = H0, limits 0..100) -> H4
%                             (-) <- [alpha] <- sqrt(max(H4,0))
%  FB:       e  = r4 - H4 -> gain Se -> Sat [-1 1] -> FIS in1
%            H4 -> (H_k - H_{k-1}) -> gain -Sde/Ts -> Sat [-1 1] -> FIS in2
%            FIS out -> gain Kdu -> accumulator I (clamp anti-windup)
%  FF:       r4 -> sqrt -> gain cFF -> F_ff
%  F1_sp  =  Sat( F_ff + I , 0, 100 )  -> setpoint of the FIC1 loop
%  ICs:      H0 = 20, F1 = cFF*sqrt(20), I = 0

function T = settleTime(y, t, ref, band, t0)
    k = find(abs(y - ref) > band, 1, 'last');
    if isempty(k),        T = 0;
    elseif k == numel(y), T = NaN;
    else,                 T = t(k+1) - t0;
    end
end