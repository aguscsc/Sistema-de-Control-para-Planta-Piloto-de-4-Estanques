%% build_FLC_L4_single.m
% LIC4 (upper tank 4) - ONE fuzzy controller (FF + P) + external integrator
%
%   Reference model (noise-free, no measured derivatives):
%     H_ref(k)  = q*H_ref(k-1) + (1-q)*r(k),      q = exp(-TsL/lambda_ref)
%     dH_ref    = (r - H_ref)/lambda_ref           slope of the trajectory
%
%   Single FLC (3 inputs, trained with ANFIS):
%     F_fis = FIS(e, H_ref, dH_ref) ~ cFF*sqrt(H_ref) + dH_ref/beta + Kp*e
%             \_____ static FF ____/   \_ overdrive _/  \_ P action _/
%     e = H_ref - H                    (tracks the trajectory, not the step)
%
%   External integrator (the only state the FIS can't hold):
%     I(k)  = I(k-1) + Ki*TsL*e        (conditional integration = anti-windup)
%     F_sp  = sat(F_fis + I, 0, 100)   -> setpoint of the FIC1 flow loop
%
%   Level loop runs at TsL = 1 s, flow loop at 0.1 s.
%
%   The closed-loop test compares:
%     1) single FLC            (this design)
%     2) ideal law             (same formula, exact - shows what the fuzzy fit costs)
%     3) static FF + I-P       (your current Simulink structure)

clear; clc; close all;

%% ================= USER SETTINGS =================
nMF        = [2 2 2];    % triangles for [e, H_ref, dH_ref]
outType    = 'constant';   % 'linear' (1st-order Sugeno) or 'constant' (singletons)
lambda_ref = 40;         % [s] reference trajectory time constant (actuator-limited, >= ~35 s)
lambda_fb  = 20;         % [s] feedback closed-loop time constant
Hd         = 30;         % [%] design point for Ti
noiseStd   = 0;          % [%] level measurement noise in the test (try 0.05)
trainPremises = false;   % false: keep uniform triangles fixed, ANFIS only fits the consequents (LSE)

%% ================= Plant (tank 4) =================
g = 9.81; A4 = 0.159; a4 = 4.1820e-04; gamma1 = 0.455;
q_max = 0.0017; h_max = 0.6; k1 = q_max;

T.beta  = (1-gamma1)*k1/(A4*h_max);                % (%/s) per % flow
T.alpha = 100*a4*sqrt(2*g*h_max/100)/(A4*h_max);   % %/s per sqrt(%)
T.cFF   = T.alpha/T.beta;                          % F_ss = cFF*sqrt(H)
T.Hmax  = (100/T.cFF)^2;                           % max steady level
T.tau_d = 2*sqrt(Hd)/T.alpha;

%% ================= Controller parameters =================
C.TsL  = 1.0;                         % [s] level loop sample time
C.lref = lambda_ref;
C.q    = exp(-C.TsL/C.lref);

C.Kp = 1/(T.beta*lambda_fb);          % % flow per % level (same at every operating point)
C.Ti = T.tau_d;                       % [s]
C.Ki = C.Kp/C.Ti;                     % % flow per (% level * s)

% Universes
C.Emax  = 10;                         % [% level] |e| beyond this: P action saturates
C.Hop   = 40;                         % [% level] H_ref in [0, Hop]   (tank max ~41.7 %)
C.DHmax = 0.4;                        % [%/s]     |dH_ref| range
C.Se    = 1/C.Emax;  C.Sh = 1/sqrt(C.Hop);  C.Sdh = 1/C.DHmax;   % Sh now scales sqrt(H)

fprintf('===== LIC4 single FLC =====\n');
fprintf('beta = %.5f  alpha = %.4f  cFF = %.3f  Hmax = %.1f %%\n', T.beta, T.alpha, T.cFF, T.Hmax);
fprintf('Kp = %.3f %%flow/%%level   Ti = %.1f s   Ki = %.4f\n', C.Kp, C.Ti, C.Ki);
fprintf('Reference: lambda_ref = %.0f s  (q = %.5f)\n', C.lref, C.q);
fprintf('Scaling:  Se = %.4f   Sh = %.4f   Sdh = %.3f\n', C.Se, C.Sh, C.Sdh);
fprintf('Max overdrive asked at a step of 10 %%: %.1f %% flow\n', 10/C.lref/T.beta);

%% ================= Training data =================
% Target in physical units (% flow); inputs normalized
targetF = @(en, hn, dhn) T.cFF*sqrt(C.Hop)*hn + (dhn*C.DHmax)/T.beta + C.Kp*(en*C.Emax);

ge = linspace(-1,1,21);  gh = linspace(0,1,21);  gd = linspace(-1,1,21);
[E3, H3, D3] = ndgrid(ge, gh, gd);
Xtr = [E3(:) H3(:) D3(:)];
Ytr = targetF(Xtr(:,1), Xtr(:,2), Xtr(:,3));

rng(1);
Xva = [2*rand(3000,1)-1, rand(3000,1), 2*rand(3000,1)-1];
Yva = targetF(Xva(:,1), Xva(:,2), Xva(:,3));

%% ================= ANFIS =================
gopt = genfisOptions('GridPartition');
gopt.NumMembershipFunctions       = nMF;
gopt.InputMembershipFunctionType  = 'trimf';
gopt.OutputMembershipFunctionType = outType;
fis0 = genfis(Xtr, Ytr, gopt);

if trainPremises
    nEp = 150;  step0 = 0.01;
else
    nEp = 2;    step0 = 1e-12;     % hybrid learning: LSE fits consequents, premises effectively frozen
end
aopt = anfisOptions('InitialFIS', fis0, 'EpochNumber', nEp, 'InitialStepSize', step0, ...
    'ValidationData', [Xva Yva], 'DisplayANFISInformation', 0, ...
    'DisplayErrorValues', 0, 'DisplayStepSize', 0, 'DisplayFinalResults', 1);
[fisT, trErr, ~, fisV, vaErr] = anfis([Xtr Ytr], aopt);

% 'warning' on purpose: gaps between trained MFs (no rule fires) must be visible
opt  = evalfisOptions('OutOfRangeInputValueMessage','none','NoRuleFiredMessage','warning');
rmseTr = @(f) sqrt(mean((evalfis(f, Xtr, opt) - Ytr).^2));
rmseVa = @(f) sqrt(mean((evalfis(f, Xva, opt) - Yva).^2));
cand = {fis0, fisT, fisV};  cname = ["initial (genfis)" "final epoch" "min validation"];
rT = cellfun(rmseTr, cand);  rV = cellfun(rmseVa, cand);
fprintf('--- RMSE of each candidate [%% flow] ---\n');
fprintf('%-18s %10s %10s\n', '', 'train', 'validation');
for i = 1:3, fprintf('%-18s %10.3f %10.3f\n', cname(i), rT(i), rV(i)); end
[~, ib] = min(rV);  fis = cand{ib};          % select on VALIDATION data (avoids overfit to the grid)
fprintf('Using: %s\n', cname(ib));

fis.Name = 'FLC_L4_single';
fis.Inputs(1).Name = 'e';  fis.Inputs(2).Name = 'Href';  fis.Inputs(3).Name = 'dHref';
fis.Outputs(1).Name = 'F';
fprintf('Rules: %d\n', numel(fis.Rules));

%% ================= Checks: slices of the surface =================
hh = linspace(0,1,101)';
Hh = (hh*sqrt(C.Hop)).^2;                 % level in % for each hn
h30 = sqrt(30/C.Hop);                     % hn that corresponds to H = 30 %
F_h  = evalfis(fis, [zeros(101,1) hh zeros(101,1)], opt);            % pure static FF
F_e  = evalfis(fis, [linspace(-1,1,101)' h30*ones(101,1) zeros(101,1)], opt);
F_dh = evalfis(fis, [zeros(101,1) h30*ones(101,1) linspace(-1,1,101)'], opt);

figure('Name','Single FLC slices','Position',[100 100 1200 360]);
subplot(1,3,1); plot(Hh, F_h, Hh, T.cFF*sqrt(Hh), '--', 'LineWidth', 1.2);
xlabel('H_{ref} [%]'); ylabel('F [%]'); title('e = 0, dH_{ref} = 0  (static FF)'); legend('FIS','cFF \surd H','Location','southeast'); grid on;
subplot(1,3,2); plot(linspace(-1,1,101)*C.Emax, F_e, 'LineWidth', 1.2); hold on;
plot(linspace(-1,1,101)*C.Emax, T.cFF*sqrt(30) + C.Kp*linspace(-1,1,101)*C.Emax, '--');
xlabel('e [%]'); title('H_{ref} = 30 %, dH_{ref} = 0  (P action)'); grid on;
subplot(1,3,3); plot(linspace(-1,1,101)*C.DHmax, F_dh, 'LineWidth', 1.2); hold on;
plot(linspace(-1,1,101)*C.DHmax, T.cFF*sqrt(30) + linspace(-1,1,101)*C.DHmax/T.beta, '--');
xlabel('dH_{ref} [%/s]'); title('H_{ref} = 30 %, e = 0  (overdrive)'); grid on;

eFF   = abs(F_h - T.cFF*sqrt(Hh));
opRng = Hh >= 10;
fprintf('Max static-FF error over 0..%d %% : %.2f %% flow\n', C.Hop, max(eFF));
fprintf('Max static-FF error over 10..%d %%: %.2f %% flow   (operating range)\n', C.Hop, max(eFF(opRng)));


% Fine-grid check over the whole input space (between training points)
[Ef, Hf, Df] = ndgrid(linspace(-1,1,41), linspace(0,1,81), linspace(-1,1,41));
Xf = [Ef(:) Hf(:) Df(:)];
errF = abs(evalfis(fis, Xf, opt) - targetF(Xf(:,1), Xf(:,2), Xf(:,3)));
fprintf('Fine-grid error: RMSE %.3f   max %.2f %% flow\n', sqrt(mean(errF.^2)), max(errF));
[~, iw] = max(errF);
[~, iw] = max(errF);
fprintf('  worst point: e_n = %.2f  H_ref = %.1f %%  dH_n = %.2f\n', Xf(iw,1), (Xf(iw,2)*sqrt(C.Hop))^2, Xf(iw,3));

% Spot checks (e = 0, dH_ref = 0 -> must equal cFF*sqrt(H))
for Hc = [10 20 30 40]
    fprintf('  F(e=0, H=%2d, dH=0) = %6.2f   expected %6.2f\n', Hc, ...
                evalfis(fis, [0 sqrt(Hc/C.Hop) 0], opt), T.cFF*sqrt(Hc));
end

%% ================= Closed-loop test =================
% Tank: Euler at 0.1 s. Flow loop approx: command ramps <= 1 %/s, then pump lag.
Ts = 0.1;  nL = round(C.TsL/Ts);
tauP = 1.3776;  aF = exp(-Ts/tauP);  rate = 1.0;
Tsim = 900;  N = round(Tsim/Ts);  t = (0:N-1)*Ts;

r  = 20*ones(1,N);  r(t >= 20) = 30;  r(t >= 550) = 15;   % up-step, then down-step
dF = zeros(1,N);    dF(t >= 300) = -5;                     % -5 % flow disturbance

rng(2);  noise = noiseStd*randn(1,N);

cases = ["single FLC", "ideal law", "static FF + I-P"];
Hs = zeros(3,N);  Fsps = zeros(3,N);
H0 = 20;  F0 = T.cFF*sqrt(H0);

for cs = 1:3
    H = H0;  F = F0;  Fc = F0;  Fsp = F0;
    Hr = H0;  I = 0;
    if cs == 3, I = C.Kp*H0; end          % I-P positional: absorbs -Kp*H at start
    for k = 1:N
        % ----- level controller (every TsL) -----
        if mod(k-1, nL) == 0
            Hm = H + noise(k);
            switch cs
                case {1,2}
                    Hr  = C.q*Hr + (1-C.q)*r(k);
                    dHr = (r(k) - Hr)/C.lref;
                    e   = Hr - Hm;
                    en  = max(-1,min(1,C.Se*e)); hn = max(0,min(1,C.Sh*sqrt(max(Hr,0))));  dn = max(-1,min(1,C.Sdh*dHr));
                    if cs == 1
                        Ffis = evalfis(fis, [en hn dn], opt);
                    else
                        Ffis = targetF(en, hn, dn);
                    end
                    base = Ffis;
                case 3
                    e    = r(k) - Hm;
                    base = T.cFF*sqrt(max(r(k),0)) - C.Kp*Hm;     % static FF + P on measurement
            end
            dI = C.Ki*C.TsL*e;
            Ftry = base + I + dI;
            if ~((Ftry > 100 && dI > 0) || (Ftry < 0 && dI < 0))
                I = I + dI;                                       % conditional integration
            end
            Fsp = min(max(base + I, 0), 100);
        end
        % ----- flow loop (approx.) -----
        Fc = Fc + max(-rate*Ts, min(rate*Ts, Fsp - Fc));
        F  = aF*F + (1-aF)*Fc;
        % ----- tank -----
        H = H + Ts*(T.beta*max(F + dF(k),0) - T.alpha*sqrt(max(H,0)));
        H = min(max(H,0),100);
        Hs(cs,k) = H;  Fsps(cs,k) = Fsp;
    end
end

figure('Name','LIC4 single FLC - closed loop');
subplot(2,1,1); plot(t, r, 'k--'); hold on; plot(t, Hs, 'LineWidth', 1.2);
legend(['setpoint' cases], 'Location','northeast'); ylabel('H_4 [%]'); grid on;
title('20\rightarrow30 % at 20 s, -5 % flow disturbance at 300 s, 30\rightarrow15 % at 550 s');
subplot(2,1,2); plot(t, Fsps, 'LineWidth', 1.1); legend(cases); ylabel('F_{sp} [%]'); xlabel('t [s]'); grid on;

fprintf('--- Closed loop (settling: 0.2 %% band) ---\n');
fprintf('%-16s %10s %10s %12s %10s %10s\n', '', 'up: OS[%]', 'up: ts[s]', 'dist: dev[%]', 'down: ts', 'Fsp std');
for cs = 1:3
    up   = t >= 20  & t < 300;
    dist = t >= 300 & t < 550;
    dn   = t >= 550;
    ts_up = settleTime(Hs(cs,up), t(up), 30, 0.2, 20);
    ts_dn = settleTime(Hs(cs,dn), t(dn), 15, 0.2, 550);
    dev   = max(abs(Hs(cs,dist) - 30));
    sd    = std(diff(Fsps(cs, t >= 250 & t < 300)));
    fprintf('%-16s %10.2f %10.0f %12.2f %10.0f %10.3f\n', cases(cs), max(Hs(cs,up)) - 30, ts_up, dev, ts_dn, sd);
end

%% ================= Save =================
writeFIS(fis, 'FLC_L4_single');
FLC_L4_single = fis;  L4S = C;  TK4 = T;
save('FLC_L4_single_params.mat', 'FLC_L4_single', 'L4S', 'TK4');

% ---------------- Simulink (level part runs at TsL = 1 s) ----------------
% r  -> Discrete Transfer Fcn  num = 1-q, den = [1 -q]  -> H_ref
% dH_ref = (r - H_ref) * (1/lambda_ref)
% e      = H_ref - H4
% e      -> gain Se  -> Sat [-1 1] -\
% H_ref  -> gain Sh  -> Sat [ 0 1] ---> Mux -> FIS (sample time 1) -> F_fis
% dH_ref -> gain Sdh -> Sat [-1 1] -/
% e -> gain Ki*TsL -> (+) -> Saturation Dynamic [-F_fis, 100-F_fis] -> I --(z^-1 back to +)
% F_sp = Sat(F_fis + I, 0, 100) -> Rate Transition (1 s -> 0.1 s) -> FIC1 setpoint
%
% ---------------- PLC (Structured Text, every 1 s) ----------------
% H_ref  := q*H_ref + (1.0-q)*r;
% dH_ref := (r - H_ref)/lambda_ref;
% e      := H_ref - H4;
% F_fis  := <FIS evaluation of LIMIT(-1,Se*e,1), LIMIT(0,Sh*H_ref,1), LIMIT(-1,Sdh*dH_ref,1)>;
% dI     := Ki*TsL*e;
% IF NOT ((F_fis+I+dI > 100.0 AND dI > 0.0) OR (F_fis+I+dI < 0.0 AND dI < 0.0)) THEN I := I + dI; END_IF;
% F1_sp  := LIMIT(0.0, F_fis + I, 100.0);

function T = settleTime(y, t, ref, band, t0)
    k = find(abs(y - ref) > band, 1, 'last');
    if isempty(k),        T = 0;
    elseif k == numel(y), T = NaN;
    else,                 T = t(k+1) - t0;
    end
end