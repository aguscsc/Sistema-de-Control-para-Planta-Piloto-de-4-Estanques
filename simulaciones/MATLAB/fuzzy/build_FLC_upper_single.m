%% build_FLC_upper_single.m
% Upper tanks LIC3 / LIC4 - ONE fuzzy controller (FF + P) + external integrator
%
%   Tank (in % units: H = level % of h_max, F = pump flow % of q_max):
%     dH/dt = beta*F - alpha*sqrt(H)
%     Tank 4: fed by pump 1 through (1-gamma1)  -> cascade on FIC1
%     Tank 3: fed by pump 2 through (1-gamma2)  -> cascade on FIC2
%
%   Reference model (noise-free, no measured derivatives):
%     H_ref(k)  = q*H_ref(k-1) + (1-q)*r(k),      q = exp(-TsL/lambda_ref)
%     dH_ref    = (r - H_ref)/lambda_ref           slope of the trajectory
%
%   Single FLC (3 inputs, trained with ANFIS, premises fixed):
%     F_fis = FIS(e, sqrt(H_ref), dH_ref) = cFF*sqrt(H_ref) + dH_ref/beta + Kp*e
%             \_______ static FF ________/   \_ overdrive _/  \_ P action _/
%     e = H_ref - H
%
%   External integrator:
%     I(k)  = I(k-1) + Ki*TsL*e        (conditional integration = anti-windup)
%     F_sp  = sat(F_fis + I, 0, 100)   -> setpoint of the flow loop
%
%   Level loop runs at TsL = 1 s, flow loop at 0.1 s.

clear; clc; close all;

%% ================= USER SETTINGS =================
TANK       = 3;          % 3 or 4
nMF        = [2 2 2];    % triangles for [e, sqrt(H_ref), dH_ref]  (2 is exact for the ideal model)
outType    = 'constant'; % singletons
lambda_ref = 40;         % [s] reference trajectory time constant (actuator-limited)
lambda_fb  = 20;         % [s] feedback closed-loop time constant
Hd         = 30;         % [%] design point for Ti
noiseStd   = 0;          % [%] level measurement noise in the test (try 0.05)
trainPremises = false;   % false: uniform triangles fixed, ANFIS only fits the consequents (LSE)
tauPump    = 1.3776;     % [s] flow-loop lag in the test (pump 2 value unknown -> pump 1 used)

%% ================= Plant =================
g = 9.81; A = 0.159; q_max = 0.0017; h_max = 0.6;
switch TANK
    case 4
        a = 4.1820e-04;  gam = 0.455;     % a4, gamma1 (pump 1)
        Hop = 40;  DHmax = 0.4;           % H_ref range [%], |dH_ref| range [%/s]
        H0 = 20;  rUp = 30;  rDn = 15;    % test: start, up-step, down-step
        loopName = 'FIC1';
    case 3
        a = 3.5752e-04;  gam = 0.4195;    % a3, gamma2 (pump 2)
        Hop = 60;  DHmax = 0.6;
        H0 = 30;  rUp = 40;  rDn = 20;
        loopName = 'FIC2';
    otherwise
        error('TANK must be 3 or 4');
end

T.beta  = (1-gam)*q_max/(A*h_max);                 % (%/s) per % flow
T.alpha = 100*a*sqrt(2*g*h_max/100)/(A*h_max);     % %/s per sqrt(%)
T.cFF   = T.alpha/T.beta;                          % F_ss = cFF*sqrt(H)
T.Hmax  = (100/T.cFF)^2;                           % max steady level
T.tau_d = 2*sqrt(Hd)/T.alpha;

%% ================= Controller parameters =================
C.TsL  = 1.0;
C.lref = lambda_ref;
C.q    = exp(-C.TsL/C.lref);

C.Kp = 1/(T.beta*lambda_fb);          % % flow per % level
C.Ti = T.tau_d;                       % [s]
C.Ki = C.Kp/C.Ti;                     % % flow per (% level * s)

C.Emax  = 10;                         % [% level]
C.Hop   = Hop;                        % [% level]  H_ref in [0, Hop]
C.DHmax = DHmax;                      % [%/s]
C.Se    = 1/C.Emax;  C.Sh = 1/sqrt(C.Hop);  C.Sdh = 1/C.DHmax;

name = sprintf('FLC_L%d_single', TANK);
fprintf('===== %s (tank %d, cascade on %s) =====\n', name, TANK, loopName);
fprintf('beta = %.5f  alpha = %.4f  cFF = %.3f  Hmax = %.1f %%\n', T.beta, T.alpha, T.cFF, T.Hmax);
fprintf('Kp = %.3f %%flow/%%level   Ti = %.1f s   Ki*TsL = %.4f\n', C.Kp, C.Ti, C.Ki*C.TsL);
fprintf('Reference: lambda_ref = %.0f s  (q = %.5f)\n', C.lref, C.q);
fprintf('Simulink gains:  e: %.4f   sqrt(H_ref): %.4f   dH_ref: %.4f (= Sdh/lambda_ref on r-H_ref)\n', ...
        C.Se, C.Sh, C.Sdh/C.lref);
fprintf('Overdrive asked at a 10 %% step: %.1f %% flow\n', 10/C.lref/T.beta);

%% ================= Training data =================
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
    nEp = 2;    step0 = 1e-12;
end
aopt = anfisOptions('InitialFIS', fis0, 'EpochNumber', nEp, 'InitialStepSize', step0, ...
    'ValidationData', [Xva Yva], 'DisplayANFISInformation', 0, ...
    'DisplayErrorValues', 0, 'DisplayStepSize', 0, 'DisplayFinalResults', 1);
[fisT, trErr, ~, fisV, vaErr] = anfis([Xtr Ytr], aopt);

opt  = evalfisOptions('OutOfRangeInputValueMessage','none','NoRuleFiredMessage','warning');
rmseTr = @(f) sqrt(mean((evalfis(f, Xtr, opt) - Ytr).^2));
rmseVa = @(f) sqrt(mean((evalfis(f, Xva, opt) - Yva).^2));
cand = {fis0, fisT, fisV};  cname = ["initial (genfis)" "final epoch" "min validation"];
rT = cellfun(rmseTr, cand);  rV = cellfun(rmseVa, cand);
fprintf('--- RMSE of each candidate [%% flow] ---\n');
fprintf('%-18s %10s %10s\n', '', 'train', 'validation');
for i = 1:3, fprintf('%-18s %10.3f %10.3f\n', cname(i), rT(i), rV(i)); end
[~, ib] = min(rV);  fis = cand{ib};
fprintf('Using: %s\n', cname(ib));

fis.Name = name;
fis.Inputs(1).Name = 'e';  fis.Inputs(2).Name = 'sqrtHref';  fis.Inputs(3).Name = 'dHref';
fis.Outputs(1).Name = 'F';
fprintf('Rules: %d\n', numel(fis.Rules));

%% ================= Checks =================
hh  = linspace(0,1,101)';
Hh  = (hh*sqrt(C.Hop)).^2;
Hm_ = 0.5*C.Hop;  hmid = sqrt(Hm_/C.Hop);        % slice point: half the range
F_h  = evalfis(fis, [zeros(101,1) hh zeros(101,1)], opt);
F_e  = evalfis(fis, [linspace(-1,1,101)' hmid*ones(101,1) zeros(101,1)], opt);
F_dh = evalfis(fis, [zeros(101,1) hmid*ones(101,1) linspace(-1,1,101)'], opt);

figure('Name',[name ' slices'],'Position',[100 100 1200 360]);
subplot(1,3,1); plot(Hh, F_h, Hh, T.cFF*sqrt(Hh), '--', 'LineWidth', 1.2);
xlabel('H_{ref} [%]'); ylabel('F [%]'); title('e = 0, dH_{ref} = 0  (static FF)');
legend('FIS','cFF \surd H','Location','southeast'); grid on;
subplot(1,3,2); plot(linspace(-1,1,101)*C.Emax, F_e, 'LineWidth', 1.2); hold on;
plot(linspace(-1,1,101)*C.Emax, T.cFF*sqrt(Hm_) + C.Kp*linspace(-1,1,101)*C.Emax, '--');
xlabel('e [%]'); title(sprintf('H_{ref} = %g %%  (P action)', Hm_)); grid on;
subplot(1,3,3); plot(linspace(-1,1,101)*C.DHmax, F_dh, 'LineWidth', 1.2); hold on;
plot(linspace(-1,1,101)*C.DHmax, T.cFF*sqrt(Hm_) + linspace(-1,1,101)*C.DHmax/T.beta, '--');
xlabel('dH_{ref} [%/s]'); title(sprintf('H_{ref} = %g %%  (overdrive)', Hm_)); grid on;

eFF = abs(F_h - T.cFF*sqrt(Hh));
fprintf('Max static-FF error over 0..%d %% : %.3f %% flow\n', C.Hop, max(eFF));

[Ef, Hf, Df] = ndgrid(linspace(-1,1,41), linspace(0,1,81), linspace(-1,1,41));
Xf = [Ef(:) Hf(:) Df(:)];
errF = abs(evalfis(fis, Xf, opt) - targetF(Xf(:,1), Xf(:,2), Xf(:,3)));
[~, iw] = max(errF);
fprintf('Fine-grid error: RMSE %.3f   max %.3f %% flow\n', sqrt(mean(errF.^2)), max(errF));
fprintf('  worst point: e_n = %.2f  H_ref = %.1f %%  dH_n = %.2f\n', Xf(iw,1), (Xf(iw,2)*sqrt(C.Hop))^2, Xf(iw,3));

for Hc = round(linspace(10, C.Hop, 4))
    fprintf('  F(e=0, H=%2d, dH=0) = %6.2f   expected %6.2f\n', Hc, ...
        evalfis(fis, [0 sqrt(Hc/C.Hop) 0], opt), T.cFF*sqrt(Hc));
end

%% ================= Closed-loop test =================
Ts = 0.1;  nL = round(C.TsL/Ts);
aF = exp(-Ts/tauPump);  rate = 1.0;
Tsim = 900;  N = round(Tsim/Ts);  t = (0:N-1)*Ts;

r  = H0*ones(1,N);  r(t >= 20) = rUp;  r(t >= 550) = rDn;
dF = zeros(1,N);    dF(t >= 300) = -5;

rng(2);  noise = noiseStd*randn(1,N);

cases = ["single FLC", "ideal law", "static FF + I-P"];
Hs = zeros(3,N);  Fsps = zeros(3,N);
F0 = T.cFF*sqrt(H0);

for cs = 1:3
    H = H0;  F = F0;  Fc = F0;  Fsp = F0;
    Hr = H0;  I = 0;
    if cs == 3, I = C.Kp*H0; end
    for k = 1:N
        if mod(k-1, nL) == 0
            Hm = H + noise(k);
            switch cs
                case {1,2}
                    Hr  = C.q*Hr + (1-C.q)*r(k);
                    dHr = (r(k) - Hr)/C.lref;
                    e   = Hr - Hm;
                    en  = max(-1,min(1,C.Se*e));
                    hn  = max(0,min(1,C.Sh*sqrt(max(Hr,0))));
                    dn  = max(-1,min(1,C.Sdh*dHr));
                    if cs == 1, base = evalfis(fis, [en hn dn], opt);
                    else,       base = targetF(en, hn, dn);
                    end
                case 3
                    e    = r(k) - Hm;
                    base = T.cFF*sqrt(max(r(k),0)) - C.Kp*Hm;
            end
            dI = C.Ki*C.TsL*e;
            Ftry = base + I + dI;
            if ~((Ftry > 100 && dI > 0) || (Ftry < 0 && dI < 0))
                I = I + dI;
            end
            Fsp = min(max(base + I, 0), 100);
        end
        Fc = Fc + max(-rate*Ts, min(rate*Ts, Fsp - Fc));
        F  = aF*F + (1-aF)*Fc;
        H = H + Ts*(T.beta*max(F + dF(k),0) - T.alpha*sqrt(max(H,0)));
        H = min(max(H,0),100);
        Hs(cs,k) = H;  Fsps(cs,k) = Fsp;
    end
end

figure('Name',[name ' - closed loop']);
subplot(2,1,1); plot(t, r, 'k--'); hold on; plot(t, Hs, 'LineWidth', 1.2);
legend(['setpoint' cases], 'Location','northeast'); ylabel(sprintf('H_%d [%%]', TANK)); grid on;
title(sprintf('%g\\rightarrow%g %% at 20 s, -5 %% flow at 300 s, %g\\rightarrow%g %% at 550 s', H0, rUp, rUp, rDn));
subplot(2,1,2); plot(t, Fsps, 'LineWidth', 1.1); legend(cases); ylabel('F_{sp} [%]'); xlabel('t [s]'); grid on;

fprintf('--- Closed loop (settling: 0.2 %% band) ---\n');
fprintf('%-16s %10s %10s %12s %10s %10s\n', '', 'up: OS[%]', 'up: ts[s]', 'dist: dev[%]', 'down: ts', 'max Fsp');
for cs = 1:3
    up   = t >= 20  & t < 300;
    dist = t >= 300 & t < 550;
    dn   = t >= 550;
    fprintf('%-16s %10.2f %10.0f %12.2f %10.0f %10.1f\n', cases(cs), ...
        max(Hs(cs,up)) - rUp, settleTime(Hs(cs,up), t(up), rUp, 0.2, 20), ...
        max(abs(Hs(cs,dist) - rUp)), settleTime(Hs(cs,dn), t(dn), rDn, 0.2, 550), max(Fsps(cs,:)));
end

%% ================= Save =================
writeFIS(fis, name);
assignin('base', name, fis);
eval(sprintf('L%dS = C; TK%d = T;', TANK, TANK));
save([name '_params.mat'], name, sprintf('L%dS', TANK), sprintf('TK%d', TANK));

% ---------------- Simulink (level part at TsL = 1 s) ----------------
% r      -> Discrete Transfer Fcn  num = 1-q, den = [1 -q]  -> H_ref
% e      =  H_ref - H                   -> gain Se            -> Sat [-1 1] -> Mux in 1
% H_ref  -> Sqrt                        -> gain Sh            -> Sat [ 0 1] -> Mux in 2
% r - H_ref                             -> gain Sdh/lambda_ref -> Sat [-1 1] -> Mux in 3
% Mux -> FIS (sample time 1) -> F_fis
% e -> gain Ki*TsL -> (+) -> Sat [-30 30] -> I  (z^-1, IC 0, back to +)
% F_sp = Sat(F_fis + I, 0, 100) -> Rate Transition (1 s -> 0.1 s) -> FIC2 (tank 3) / FIC1 (tank 4)

function T = settleTime(y, t, ref, band, t0)
    k = find(abs(y - ref) > band, 1, 'last');
    if isempty(k),        T = 0;
    elseif k == numel(y), T = NaN;
    else,                 T = t(k+1) - t0;
    end
end