%% build_FLC_flow_anfis.m
% Flow loops FIC1 / FIC2 - fuzzy controller trained with ANFIS
%
% Target control law (velocity form, keeps integral action):
%   s  = e + Ti*de                 predicted error if the CV were frozen now
%   du = Kdu * sat(s / delta)      [% per sample]
%   - |s| < delta : exactly the Dahlin PI
%   - |s| > delta : full rate (pump rate limit), braking as s -> delta
%
% Normalized inputs : e_n = Se*e, de_n = Sde*de   (saturated to [-1,1] in Simulink)
% Normalized output : du_n in [-1,1],  du = Kdu*du_n
% de ON MEASUREMENT : de = -(y_k - y_{k-1})/Ts
%
% FIS: zero-order Sugeno (singletons), 3x3 triangular MFs
%      optional reduction to nOut shared singletons (NB..PB)

clear; clc; close all;

%% ================= USER SETTINGS =================
PUMP   = 1;                      % 1 or 2
nMF    = 3;                      % triangles per input
nOut   = 5;                      % shared singletons after training (0 = keep one per rule)

tauP   = [1.3776, NaN];          % slowest time constant of each pump [s]  <-- fill pump 2

%% ================= Parameters =================
P.Ts     = 0.1;                  % [s]   PLC sample time
P.K      = 1.0;                  % plant gain through the PLC mapping
P.tau    = tauP(PUMP);
P.rate   = 1.0;                  % [%/s] pump rate limit
if isnan(P.tau), error('Fill in tauP(%d) first.', PUMP); end
P.lambda = P.tau/2;              % [s]   Dahlin closed-loop time constant

% Dahlin PI
P.a  = exp(-P.Ts/P.tau);
P.q  = exp(-P.Ts/P.lambda);
P.Kc = (1 - P.q)/(P.K*(1 - P.a));
P.Ti = P.Ts/(1 - P.a);

% Control law
P.Kdu   = P.rate*P.Ts;                     % [%/sample]
P.delta = P.Kdu*P.Ti/(P.Kc*P.Ts);          % [%] linear band (matches Dahlin PI)

% Universes and scaling
P.DEmax = 1.25*P.rate;                     % [%/s]
P.Emax  = P.Ti*P.DEmax + 2*P.delta;        % [%]
P.Se    = 1/P.Emax;
P.Sde   = 1/P.DEmax;

name = sprintf('FLC_F%d', PUMP);
fprintf('===== %s =====\n', name);
fprintf('--- Dahlin PI ---\n');
fprintf('Kc = %.4f   Ti = %.4f s\n', P.Kc, P.Ti);
fprintf('--- Control law ---\n');
fprintf('Kdu = %.3f %%/sample   delta = %.3f %%\n', P.Kdu, P.delta);
fprintf('--- Universes / scaling ---\n');
fprintf('Emax  = %.3f %%    -> Se  = %.4f\n', P.Emax,  P.Se);
fprintf('DEmax = %.3f %%/s  -> Sde = %.4f   (Simulink gain on -dy: Sde/Ts = %.3f)\n', ...
        P.DEmax, P.Sde, P.Sde/P.Ts);

%% ================= Training data =================
target = @(en, den) max(-1, min(1, (en*P.Emax + P.Ti*den*P.DEmax)/P.delta));

g = linspace(-1, 1, 61);
[EN, DEN] = meshgrid(g);
Xtr = [EN(:) DEN(:)];
Ytr = target(Xtr(:,1), Xtr(:,2));

rng(1);
Xva = 2*rand(1500,2) - 1;
Yva = target(Xva(:,1), Xva(:,2));

%% ================= ANFIS =================
gopt = genfisOptions('GridPartition');
gopt.NumMembershipFunctions       = [nMF nMF];
gopt.InputMembershipFunctionType  = 'trimf';
gopt.OutputMembershipFunctionType = 'constant';     % singletons
fis0 = genfis(Xtr, Ytr, gopt);

aopt = anfisOptions('InitialFIS', fis0, ...
    'EpochNumber', 200, ...
    'ValidationData', [Xva Yva], ...
    'DisplayANFISInformation', 0, ...
    'DisplayErrorValues', 0, ...
    'DisplayStepSize', 0, ...
    'DisplayFinalResults', 1);
[fisT, trErr, ~, fisV, vaErr] = anfis([Xtr Ytr], aopt);

% Evaluate every candidate ourselves and keep the best
opt  = evalfisOptions('OutOfRangeInputValueMessage','none','NoRuleFiredMessage','none');
rmse = @(f) sqrt(mean((evalfis(f, Xtr, opt) - Ytr).^2));
cand  = {fis0, fisT, fisV};
cname = ["initial (genfis)" "final epoch" "min validation"];
r = cellfun(rmse, cand);
fprintf('--- evalfis RMSE of each candidate ---\n');
for i = 1:3, fprintf('%-18s %.4f\n', cname(i), r(i)); end
[~, ib] = min(r);  fis = cand{ib};
fprintf('Using: %s\n', cname(ib));

fis.Name = name;
fis.Inputs(1).Name  = 'e';
fis.Inputs(2).Name  = 'de';
fis.Outputs(1).Name = 'du';

% Trained singletons as a rule table (rows: e, cols: de)
nR = numel(fis.Rules);
cTab = nan(nMF);
for k = 1:nR
    ant = fis.Rules(k).Antecedent;
    cTab(ant(1), ant(2)) = fis.Outputs(1).MembershipFunctions(fis.Rules(k).Consequent).Parameters;
end
fprintf('--- Trained singletons (rows e = N..P, cols de = N..P) ---\n');
disp(round(cTab, 3));

%% ================= Optional: reduce to nOut shared singletons =================
if nOut > 0
    c = zeros(nR,1);
    for k = 1:nR
        c(k) = fis.Outputs(1).MembershipFunctions(fis.Rules(k).Consequent).Parameters;
    end

    % 1) label each rule with the nearest of nOut levels spread over the trained values
    lev0   = linspace(min(c), max(c), nOut);
    [~, L] = min(abs(c - lev0), [], 2);

    % 2) normalized firing strengths on the training grid (premises fixed)
    W = zeros(size(Xtr,1), nR);
    for k = 1:size(Xtr,1)
        [~,~,~,~,rf] = evalfis(fis, Xtr(k,:), opt);
        W(k,:) = rf(:)' / max(sum(rf), eps);
    end

    % 3) least squares for the level values
    A = zeros(size(Xtr,1), nOut);
    for j = 1:nOut, A(:,j) = sum(W(:, L == j), 2); end
    used = any(A, 1);
    v = lev0;  v(used) = A(:,used) \ Ytr;

    % 4) rebuild a fresh Sugeno FIS
    if nOut == 5
        outNames = ["NB" "NS" "ZE" "PS" "PB"];
    elseif nOut == 3
        outNames = ["N" "Z" "P"];
    else
        outNames = "o" + string(1:nOut);
    end
    fisR = sugfis('Name', name, 'AndMethod','prod', 'DefuzzificationMethod','wtaver');
    for i = 1:2
        fisR = addInput(fisR, [-1 1], 'Name', fis.Inputs(i).Name);
        for m = 1:numel(fis.Inputs(i).MembershipFunctions)
            mf = fis.Inputs(i).MembershipFunctions(m);
            fisR = addMF(fisR, fis.Inputs(i).Name, mf.Type, mf.Parameters, 'Name', mf.Name);
        end
    end
    fisR = addOutput(fisR, [min(-1,min(v)) max(1,max(v))], 'Name', 'du');
    for j = 1:nOut
        fisR = addMF(fisR, 'du', 'constant', v(j), 'Name', outNames(j));
    end
    ruleList = zeros(nR, 5);
    for k = 1:nR
        ruleList(k,:) = [fis.Rules(k).Antecedent, L(k), 1, 1];   % [in1 in2 out weight AND]
    end
    fisR = addRule(fisR, ruleList);

    fprintf('--- Reduction 9 -> %d shared singletons ---\n', nOut);
    fprintf('LS fit RMSE: %.4f\n', sqrt(mean((A*v(:) - Ytr).^2)));
    fprintf('Singletons : %s\n', mat2str(v, 3));
    lTab = nan(nMF);
    for k = 1:nR, ant = fis.Rules(k).Antecedent; lTab(ant(1),ant(2)) = L(k); end
    fprintf('Rule table (label index, rows e = N..P, cols de = N..P):\n');
    disp(lTab);

    fis = fisR;
end

%% ================= Checks =================
DU = reshape(evalfis(fis, Xtr, opt), size(EN));
TG = reshape(Ytr, size(EN));
fprintf('--- Fit ---\n');
fprintf('RMSE train = %.4f   max |err| = %.4f\n', sqrt(mean((DU(:)-TG(:)).^2)), max(abs(DU(:)-TG(:))));

h  = 1e-3;
ge = (evalfis(fis,[h 0],opt) - evalfis(fis,[-h 0],opt))/(2*h);
gd = (evalfis(fis,[0 h],opt) - evalfis(fis,[0 -h],opt))/(2*h);
fprintf('Origin slope e : %.2f (target %.2f)   de : %.2f (target %.2f)\n', ...
        ge, P.Emax/P.delta, gd, P.Ti*P.DEmax/P.delta);

figure('Name',[name ' surface'],'Position',[100 100 1200 380]);
subplot(1,3,1); surf(EN,DEN,TG,'EdgeColor','none'); view(-35,30);
title('Target law'); xlabel('e_n'); ylabel('de_n'); zlabel('du_n');
subplot(1,3,2); surf(EN,DEN,DU,'EdgeColor','none'); view(-35,30);
title('FIS'); xlabel('e_n'); ylabel('de_n'); zlabel('du_n');
subplot(1,3,3); contourf(EN,DEN,DU-TG,20,'LineColor','none'); colorbar;
title('Error'); xlabel('e_n'); ylabel('de_n'); axis square;

figure('Name','Training'); plot(trErr,'LineWidth',1.2); hold on; plot(vaErr,'LineWidth',1.2);
legend('train','validation'); xlabel('epoch'); ylabel('RMSE'); grid on;

figure('Name','Input MFs');
subplot(2,1,1); plotmf(fis,'input',1); title('e_n');
subplot(2,1,2); plotmf(fis,'input',2); title('de_n');

%% ================= Quick closed-loop test =================
% Plant 1/(tau s + 1) with ZOH, unity gain, one-sample PLC delay, CV in [0,100]
T  = 60;  N = round(T/P.Ts);  t = (0:N-1)*P.Ts;
r  = 40*ones(1,N);  r(t >= 10) = 50;

y_f = 40*ones(1,N); u_f = 40*ones(1,N); uo_f = 40;
y_p = 40*ones(1,N); u_p = 40*ones(1,N); uo_p = 40;
for k = 2:N
    y_f(k) = P.a*y_f(k-1) + (1-P.a)*P.K*uo_f;
    y_p(k) = P.a*y_p(k-1) + (1-P.a)*P.K*uo_p;

    % FLC
    e   = r(k) - y_f(k);
    de  = -(y_f(k) - y_f(k-1))/P.Ts;
    en  = max(-1, min(1, P.Se*e));
    den = max(-1, min(1, P.Sde*de));
    du  = P.Kdu*max(-1, min(1, evalfis(fis,[en den],opt)));
    uo_f   = u_f(k-1);
    u_f(k) = max(0, min(100, u_f(k-1) + du));

    % Dahlin PI (P on measurement), same actuator clamp
    e   = r(k) - y_p(k);
    du  = -P.Kc*(y_p(k) - y_p(k-1)) + P.Kc*(1-P.a)*e;
    du  = max(-P.Kdu, min(P.Kdu, du));
    uo_p   = u_p(k-1);
    u_p(k) = max(0, min(100, u_p(k-1) + du));
end

figure('Name',[name ': step 40 -> 50 %']);
subplot(2,1,1);
plot(t, r, 'k--', t, y_f, 'LineWidth', 1.3); hold on; plot(t, y_p, 'LineWidth', 1.3);
legend('setpoint','FLC','Dahlin PI + rate clamp','Location','southeast');
ylabel('flow [%]'); grid on;
subplot(2,1,2);
plot(t, u_f, t, u_p, 'LineWidth', 1.3);
legend('CV FLC','CV PI','Location','southeast'); ylabel('CV [%]'); xlabel('t [s]'); grid on;

fprintf('--- Step 40->50 ---\n');
fprintf('FLC: overshoot %.2f %%   settling %.1f s\n', max(y_f) - 50, settleTime(y_f, t, 50, 0.2, 10));
fprintf('PI : overshoot %.2f %%   settling %.1f s\n', max(y_p) - 50, settleTime(y_p, t, 50, 0.2, 10));
fprintf('Ideal (rate-limited) ~ %.1f s\n', 10/P.rate + P.tau);

%% ================= Save for Simulink =================
writeFIS(fis, name);                       % <name>.fis
assignin('base', name, fis);               % FLC_F1 / FLC_F2 in the workspace
eval(sprintf('P%d = P;', PUMP));           % P1 / P2
save([name '_params.mat'], name, sprintf('P%d', PUMP));

% Simulink:
%  e  -> gain Se            -> Sat [-1 1] -> FIS in1
%  y  -> (y_k - y_{k-1})    -> gain -Sde/Ts -> Sat [-1 1] -> FIS in2
%  FIS out -> gain Kdu (0.1) -> clamp [-0.1 0.1] -> accumulator -> Sat [0 100] -> CV
%  Accumulator and plant ICs at the operating point (e.g. 40).

%% ================= Local functions =================
function T = settleTime(y, t, ref, band, t0)
    k = find(abs(y - ref) > band, 1, 'last');
    if isempty(k),        T = 0;
    elseif k == numel(y), T = NaN;          % never settled
    else,                 T = t(k+1) - t0;
    end
end