%% build_FLC_lower_single.m
% Lower tanks LIC1 / LIC2 - ONE fuzzy controller (FF + P) + external integrator
%
%   Tank (in % units: H = level % of h_max, F = pump flow % of q_max):
%     dH/dt = beta*F - alpha*sqrt(H) + alphaU*sqrt(H_up)
%     Tank 1: pump 1 (gamma1 fraction), inflow from tank 3   -> cascade on FIC1
%     Tank 2: pump 2 (gamma2 fraction), inflow from tank 4   -> cascade on FIC2
%   (LIC1/LIC2 drive the pumps directly; the upper tank is a measured disturbance.)
%
%   Inverse model (what flow makes the level follow H_ref):
%     F = c1*sqrt(H_ref) - cU*sqrt(H_up) + dH_ref/beta
%         \_ hold level _/ \_ upper inflow_/ \_overdrive_/
%
%   Single FLC (4 inputs, trained with ANFIS, premises fixed):
%     F_fis = FIS(e, sqrt(H_ref), dH_ref, sqrt(H_up)) = c1*sqrt(H_ref) - cU*sqrt(H_up) + dH_ref/beta + Kp*e
%   External integrator:
%     I(k) = I(k-1) + Ki*TsL*e,   F_sp = sat(F_fis + I, 0, 100)
%
%   The closed-loop test compares:
%     1) single FLC
%     2) ideal law                     (exact formula)
%     3) FLC without upper-tank FF     (sqrt(H_up) frozen at its initial value)

clear; clc; close all;

%% ================= USER SETTINGS =================
TANK       = 2;            % 1 or 2
nMF        = [2 2 2 2];    % triangles for [e, sqrt(H_ref), dH_ref, sqrt(H_up)]  (2 is exact)
outType    = 'constant';
lambda_ref = 40;           % [s]
lambda_fb  = 20;           % [s]
Hd         = 20;           % [%] design point for Ti
noiseStd   = 0;            % [%] noise on both level measurements in the test
trainPremises = false;
tauPump    = 1.3776;       % [s] flow-loop lag in the test

%% ================= Plant =================
g = 9.81; A = 0.159; q_max = 0.0017; h_max = 0.6;
a1 = 6.1360e-04; a2 = 5.4678e-04; a3 = 3.5752e-04; a4 = 4.1820e-04;
gamma1 = 0.455;  gamma2 = 0.4195;

switch TANK
    case 1
        a = a1;  gam = gamma1;          % own outlet, own pump fraction (pump 1)
        aU = a3; gamU = gamma2;         % upper tank 3, fed by (1-gamma2) of pump 2
        H0 = 20;  HU0 = 30;             % test: initial levels
        loopName = 'FIC1';  upName = 'tank 3';
    case 2
        a = a2;  gam = gamma2;
        aU = a4; gamU = gamma1;         % upper tank 4, fed by (1-gamma1) of pump 1
        H0 = 20;  HU0 = 20;
        loopName = 'FIC2';  upName = 'tank 4';
    otherwise
        error('TANK must be 1 or 2');
end
rUp = 30;  rDn = 20;  dFup = -10;       % test: up-step, down-step, upper-pump change [% flow]

s = sqrt(2*g*h_max/100);
T.beta   = gam*q_max/(A*h_max);
T.alpha  = 100*a *s/(A*h_max);
T.alphaU = 100*aU*s/(A*h_max);
T.c1     = T.alpha/T.beta;               % flow to hold the level (no upper inflow)
T.cU     = T.alphaU/T.beta;              % flow equivalent of the upper-tank inflow
T.tau_d  = 2*sqrt(Hd)/T.alpha;

U.beta  = (1-gamU)*q_max/(A*h_max);      % upper tank (only for the test)
U.alpha = T.alphaU;
U.cFF   = U.alpha/U.beta;

%% ================= Controller parameters =================
C.TsL  = 1.0;
C.lref = lambda_ref;
C.q    = exp(-C.TsL/C.lref);
C.Kp = 1/(T.beta*lambda_fb);
C.Ti = T.tau_d;
C.Ki = C.Kp/C.Ti;

C.Emax  = 10;      C.Hop = 50;      C.DHmax = 0.5;     C.HUop = 60;
C.Se  = 1/C.Emax;  C.Sh  = 1/sqrt(C.Hop);  C.Sdh = 1/C.DHmax;  C.Su = 1/sqrt(C.HUop);

name = sprintf('FLC_L%d_single', TANK);
fprintf('===== %s (tank %d, pump via %s, upper inflow from %s) =====\n', name, TANK, loopName, upName);
fprintf('beta = %.5f  alpha = %.4f  alphaU = %.4f\n', T.beta, T.alpha, T.alphaU);
fprintf('Steady state:  F = %.2f*sqrt(H) - %.2f*sqrt(H_up)\n', T.c1, T.cU);
Hlo = ((T.cU*sqrt(HU0))/T.c1)^2;  Hhi = ((100 + T.cU*sqrt(HU0))/T.c1)^2;
fprintf('Reachable level with H_up = %g %%:  %.1f %% (pump 0) .. %.1f %% (pump 100)\n', HU0, Hlo, Hhi);
fprintf('Kp = %.3f   Ti = %.1f s   Ki*TsL = %.4f\n', C.Kp, C.Ti, C.Ki*C.TsL);
fprintf('Simulink gains:  e: %.4f   sqrt(H_ref): %.4f   r-H_ref: %.4f   sqrt(H_up): %.4f\n', ...
        C.Se, C.Sh, C.Sdh/C.lref, C.Su);
fprintf('Overdrive asked at a 10 %% step: %.1f %% flow\n', 10/C.lref/T.beta);

%% ================= Training data =================
targetF = @(en, hn, dn, un) T.c1*sqrt(C.Hop)*hn - T.cU*sqrt(C.HUop)*un ...
                            + (dn*C.DHmax)/T.beta + C.Kp*(en*C.Emax);
g1 = linspace(-1,1,11);  g0 = linspace(0,1,11);
[E4, H4, D4, U4] = ndgrid(g1, g0, g1, g0);
Xtr = [E4(:) H4(:) D4(:) U4(:)];
Ytr = targetF(Xtr(:,1), Xtr(:,2), Xtr(:,3), Xtr(:,4));
rng(1);
Xva = [2*rand(3000,1)-1, rand(3000,1), 2*rand(3000,1)-1, rand(3000,1)];
Yva = targetF(Xva(:,1), Xva(:,2), Xva(:,3), Xva(:,4));

%% ================= ANFIS =================
gopt = genfisOptions('GridPartition');
gopt.NumMembershipFunctions       = nMF;
gopt.InputMembershipFunctionType  = 'trimf';
gopt.OutputMembershipFunctionType = outType;
fis0 = genfis(Xtr, Ytr, gopt);
if trainPremises, nEp = 150; step0 = 0.01; else, nEp = 2; step0 = 1e-12; end
aopt = anfisOptions('InitialFIS', fis0, 'EpochNumber', nEp, 'InitialStepSize', step0, ...
    'ValidationData', [Xva Yva], 'DisplayANFISInformation', 0, ...
    'DisplayErrorValues', 0, 'DisplayStepSize', 0, 'DisplayFinalResults', 1);
[fisT, ~, ~, fisV] = anfis([Xtr Ytr], aopt);

opt  = evalfisOptions('OutOfRangeInputValueMessage','none','NoRuleFiredMessage','warning');
rmse = @(f, X, Y) sqrt(mean((evalfis(f, X, opt) - Y).^2));
cand = {fis0, fisT, fisV};  cname = ["initial (genfis)" "final epoch" "min validation"];
rT = cellfun(@(f) rmse(f, Xtr, Ytr), cand);  rV = cellfun(@(f) rmse(f, Xva, Yva), cand);
fprintf('--- RMSE of each candidate [%% flow] ---\n');
for i = 1:3, fprintf('%-18s train %8.3f   validation %8.3f\n', cname(i), rT(i), rV(i)); end
[~, ib] = min(rV);  fis = cand{ib};  fprintf('Using: %s\n', cname(ib));

fis.Name = name;
fis.Inputs(1).Name = 'e';  fis.Inputs(2).Name = 'sqrtHref';
fis.Inputs(3).Name = 'dHref';  fis.Inputs(4).Name = 'sqrtHup';
fis.Outputs(1).Name = 'F';
fprintf('Rules: %d\n', numel(fis.Rules));

%% ================= Checks =================
fprintf('--- Spot checks (e = 0, dH_ref = 0) ---\n');
for Hc = [15 20 30]
    for Hu = [10 30 50]
        Fx = T.c1*sqrt(Hc) - T.cU*sqrt(Hu);
        fprintf('  H = %2d, H_up = %2d:  FIS %7.2f   expected %7.2f\n', Hc, Hu, ...
            evalfis(fis, [0 sqrt(Hc/C.Hop) 0 sqrt(Hu/C.HUop)], opt), Fx);
    end
end

%% ================= Closed-loop test =================
Ts = 0.1;  nL = round(C.TsL/Ts);  aF = exp(-Ts/tauPump);  rate = 1.0;
Tsim = 900;  N = round(Tsim/Ts);  t = (0:N-1)*Ts;
r    = H0*ones(1,N);  r(t >= 20) = rUp;  r(t >= 550) = rDn;
FU0  = U.cFF*sqrt(HU0);
FUsp = FU0*ones(1,N);  FUsp(t >= 300) = FU0 + dFup;       % upper tank's pump changes at 300 s
rng(2);  nz = noiseStd*randn(2,N);

cases = ["single FLC", "ideal law", "no upper-tank FF"];
Hs = zeros(3,N);  HUs = zeros(3,N);  Fsps = zeros(3,N);
F0 = T.c1*sqrt(H0) - T.cU*sqrt(HU0);
un0 = C.Su*sqrt(HU0);

for cs = 1:3
    H = H0;  F = F0;  Fc = F0;  Fsp = F0;  Hr = H0;  I = 0;
    HU = HU0;  FU = FU0;  FUc = FU0;
    for k = 1:N
        if mod(k-1, nL) == 0
            Hm = H + nz(1,k);  HUm = HU + nz(2,k);
            Hr  = C.q*Hr + (1-C.q)*r(k);
            dHr = (r(k) - Hr)/C.lref;
            e   = Hr - Hm;
            en  = max(-1,min(1,C.Se*e));
            hn  = max(0,min(1,C.Sh*sqrt(max(Hr,0))));
            dn  = max(-1,min(1,C.Sdh*dHr));
            un  = max(0,min(1,C.Su*sqrt(max(HUm,0))));
            switch cs
                case 1, base = evalfis(fis, [en hn dn un], opt);
                case 2, base = targetF(en, hn, dn, un);
                case 3, base = evalfis(fis, [en hn dn un0], opt);
            end
            dI = C.Ki*C.TsL*e;
            Ftry = base + I + dI;
            if ~((Ftry > 100 && dI > 0) || (Ftry < 0 && dI < 0)), I = I + dI; end
            Fsp = min(max(base + I, 0), 100);
        end
        % own pump (flow loop approx.)
        Fc = Fc + max(-rate*Ts, min(rate*Ts, Fsp - Fc));   F = aF*F + (1-aF)*Fc;
        % upper tank's pump (other loop, here just a setpoint change)
        FUc = FUc + max(-rate*Ts, min(rate*Ts, FUsp(k) - FUc));  FU = aF*FU + (1-aF)*FUc;
        % tanks
        HU = HU + Ts*(U.beta*max(FU,0) - U.alpha*sqrt(max(HU,0)));      HU = min(max(HU,0),100);
        H  = H  + Ts*(T.beta*max(F,0) - T.alpha*sqrt(max(H,0)) + T.alphaU*sqrt(max(HU,0)));
        H  = min(max(H,0),100);
        Hs(cs,k) = H;  HUs(cs,k) = HU;  Fsps(cs,k) = Fsp;
    end
end

figure('Name',[name ' - closed loop']);
subplot(3,1,1); plot(t, r, 'k--'); hold on; plot(t, Hs, 'LineWidth', 1.2);
legend(['setpoint' cases], 'Location','northeast'); ylabel(sprintf('H_%d [%%]', TANK)); grid on;
title(sprintf('%g\\rightarrow%g %% at 20 s, upper pump %+g %% at 300 s, %g\\rightarrow%g %% at 550 s', H0, rUp, dFup, rUp, rDn));
subplot(3,1,2); plot(t, HUs(1,:), 'LineWidth', 1.2); ylabel(['H ' upName ' [%]']); grid on;
subplot(3,1,3); plot(t, Fsps, 'LineWidth', 1.1); legend(cases); ylabel('F_{sp} [%]'); xlabel('t [s]'); grid on;

fprintf('--- Closed loop (settling: 0.2 %% band) ---\n');
fprintf('%-18s %10s %10s %12s %10s %8s %8s\n', '', 'up: OS[%]', 'up: ts[s]', 'dist: dev[%]', 'down: ts', 'min Fsp', 'max Fsp');
for cs = 1:3
    up = t >= 20 & t < 300;  dist = t >= 300 & t < 550;  dn = t >= 550;
    fprintf('%-18s %10.2f %10.0f %12.2f %10.0f %8.1f %8.1f\n', cases(cs), ...
        max(Hs(cs,up)) - rUp, settleTime(Hs(cs,up), t(up), rUp, 0.2, 20), ...
        max(abs(Hs(cs,dist) - rUp)), settleTime(Hs(cs,dn), t(dn), rDn, 0.2, 550), ...
        min(Fsps(cs,:)), max(Fsps(cs,:)));
end

%% ================= Save =================
writeFIS(fis, name);
assignin('base', name, fis);
eval(sprintf('L%dS = C; TK%d = T;', TANK, TANK));
save([name '_params.mat'], name, sprintf('L%dS', TANK), sprintf('TK%d', TANK));

% ---------------- Simulink (level part at TsL = 1 s) ----------------
% Same as the upper tanks, plus a 4th FIS input:
%   r      -> Discrete Transfer Fcn (1-q)/(1-q z^-1)            -> H_ref
%   e      =  H_ref - H              -> gain Se             -> Sat [-1 1] -> Mux in 1
%   H_ref  -> Sqrt                   -> gain Sh             -> Sat [ 0 1] -> Mux in 2
%   r - H_ref                        -> gain Sdh/lambda_ref -> Sat [-1 1] -> Mux in 3
%   H_up   -> Sqrt (of max(H_up,0))  -> gain Su             -> Sat [ 0 1] -> Mux in 4
%   Mux -> FIS (sample time 1) -> F_fis;   integrator on e (gain Ki*TsL), F_sp = Sat(F_fis + I, 0, 100)
%   H_up = level of tank 3 (for LIC1) or tank 4 (for LIC2)

function T = settleTime(y, t, ref, band, t0)
    k = find(abs(y - ref) > band, 1, 'last');
    if isempty(k),        T = 0;
    elseif k == numel(y), T = NaN;
    else,                 T = t(k+1) - t0;
    end
end