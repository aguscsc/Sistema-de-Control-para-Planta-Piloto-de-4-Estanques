%% build_FLC_F1_anfis.m
% Flow loop FIC1 (pump 1) - Phase 3 via ANFIS
%
% ANFIS learns the target control law (velocity form, so integral action is kept):
%
%   s     = e + Ti*de                 predicted error if CV were frozen now
%   du    = Kdu * sat(s / delta)      [% per sample]
%
%   - Near the setpoint (|s| < delta) it is EXACTLY the Dahlin PI
%   - Far from it, it ramps at the pump rate limit (Kdu = rate*Ts)
%   - It starts braking when s -> delta, so it does not overshoot
%
% Inputs are normalized:  e_n = e/Emax,  de_n = de/DEmax  (saturated to [-1,1])
% Output is normalized:   du_n in [-1,1],  du = Kdu*du_n
% de is computed ON MEASUREMENT:  de = -(y_k - y_{k-1})/Ts
%
% The final section simulates the loop in MATLAB (FLC vs Dahlin PI)
% before you take it to Simulink.

clear; clc; close all;

%% ---------------- Parameters ----------------
P.Ts     = 0.1;          % [s]   PLC sample time
P.K      = 1.0;          % plant gain seen through the PLC mapping
P.tau    = 1.3776;       % [s]   slowest pump time constant (design case)
P.rate   = 1.0;          % [%/s] pump rate limit
P.lambda = P.tau/2;      % [s]   Dahlin closed-loop time constant near setpoint

% Dahlin PI
P.a  = exp(-P.Ts/P.tau);
P.q  = exp(-P.Ts/P.lambda);
P.Kc = (1 - P.q)/(P.K*(1 - P.a));
P.Ti = P.Ts/(1 - P.a);

% Control law
P.Kdu   = P.rate*P.Ts;                    % [%/sample] max move = rate limit
P.delta = P.Kdu*P.Ti/(P.Kc*P.Ts);         % [%] linear band -> matches Dahlin PI exactly

% Input universes (physical), chosen so saturation never hides the sign of s
P.DEmax = 1.25*P.rate;                    % [%/s] flow can't change faster than the pump
P.Emax  = P.Ti*P.DEmax + 2*P.delta;       % [%]   at (e=+Emax, de=-DEmax) s is still >= delta

% Scaling factors for Simulink
P.Se  = 1/P.Emax;                         % e_n  = Se*e
P.Sde = 1/P.DEmax;                        % de_n = Sde*de,  de in %/s

fprintf('--- Dahlin PI ---\n');
fprintf('Kc = %.4f   Ti = %.4f s\n', P.Kc, P.Ti);
fprintf('--- Control law ---\n');
fprintf('Kdu = %.3f %%/sample   delta = %.3f %%\n', P.Kdu, P.delta);
fprintf('--- Universes / scaling ---\n');
fprintf('Emax  = %.3f %%    -> Se  = %.4f\n', P.Emax,  P.Se);
fprintf('DEmax = %.3f %%/s  -> Sde = %.4f   (Simulink gain on -dy: Sde/Ts = %.3f)\n', ...
        P.DEmax, P.Sde, P.Sde/P.Ts);

%% ---------------- Training data ----------------
target = @(en, den) max(-1, min(1, (en*P.Emax + P.Ti*den*P.DEmax)/P.delta));

% Training: dense grid
g = linspace(-1, 1, 61);
[EN, DEN] = meshgrid(g);
Xtr = [EN(:) DEN(:)];
Ytr = target(Xtr(:,1), Xtr(:,2));

% Validation: random points
rng(1);
Xva = 2*rand(1500,2) - 1;
Yva = target(Xva(:,1), Xva(:,2));

%% ---------------- ANFIS ----------------
gopt = genfisOptions('GridPartition');
gopt.NumMembershipFunctions       = [3 3];      % was [5 5]
gopt.InputMembershipFunctionType  = 'trimf';    % was 'gaussmf'
gopt.OutputMembershipFunctionType = 'constant';
fis0 = genfis(Xtr, Ytr, gopt);

aopt = anfisOptions('InitialFIS', fis0, ...
    'EpochNumber', 200, ...
    'ValidationData', [Xva Yva], ...
    'DisplayANFISInformation', 0, ...
    'DisplayErrorValues', 0, ...
    'DisplayStepSize', 0, ...
    'DisplayFinalResults', 1);
[fisT, trErr, ~, fisV, vaErr] = anfis([Xtr Ytr], aopt);

% Sanity check: evaluate every candidate ourselves and keep the best one
opt  = evalfisOptions('OutOfRangeInputValueMessage','none','NoRuleFiredMessage','none');
rmse = @(f) sqrt(mean((evalfis(f, Xtr, opt) - Ytr).^2));
cand = {fis0, fisT, fisV};  names = ["initial (genfis)" "final epoch" "min validation"];
r    = cellfun(rmse, cand);
fprintf('--- evalfis RMSE of each candidate ---\n');
for i = 1:3, fprintf('%-18s %.4f\n', names(i), r(i)); end
[~, ib] = min(r);  fis = cand{ib};
fprintf('Using: %s\n', names(ib));

fis.Name = 'FLC_F1';
fis.Inputs(1).Name = 'e';
fis.Inputs(2).Name = 'de';
fis.Outputs(1).Name = 'du';
for i = 1:2
    fis.Inputs(i).Range = [-1 1];
end

%% ---------------- Checks ----------------
DU  = reshape(evalfis(fis, Xtr, opt), size(EN));
TG  = reshape(Ytr, size(EN));
fprintf('--- Fit ---\n');
fprintf('RMSE train = %.4f   max |err| = %.4f\n', sqrt(mean((DU(:)-TG(:)).^2)), max(abs(DU(:)-TG(:))));

% Gain near the origin vs Dahlin PI (should be close to 1)
h  = 1e-3;
ge = (evalfis(fis,[h 0],opt) - evalfis(fis,[-h 0],opt))/(2*h);   % d du_n / d e_n
gd = (evalfis(fis,[0 h],opt) - evalfis(fis,[0 -h],opt))/(2*h);   % d du_n / d de_n
ge_target = P.Emax/P.delta;
gd_target = P.Ti*P.DEmax/P.delta;
fprintf('Origin slope e : %.2f (target %.2f)   de : %.2f (target %.2f)\n', ge, ge_target, gd, gd_target);

figure('Name','FLC_F1 ANFIS','Position',[100 100 1200 380]);
subplot(1,3,1); surf(EN,DEN,TG,'EdgeColor','none'); view(-35,30);
title('Target law'); xlabel('e_n'); ylabel('de_n'); zlabel('du_n');
subplot(1,3,2); surf(EN,DEN,DU,'EdgeColor','none'); view(-35,30);
title('ANFIS'); xlabel('e_n'); ylabel('de_n'); zlabel('du_n');
subplot(1,3,3); contourf(EN,DEN,DU-TG,20,'LineColor','none'); colorbar;
title('Error'); xlabel('e_n'); ylabel('de_n'); axis square;

figure('Name','Training'); plot(trErr,'LineWidth',1.2); hold on; plot(vaErr,'LineWidth',1.2);
legend('train','validation'); xlabel('epoch'); ylabel('RMSE'); grid on;

figure('Name','Input MFs');
subplot(2,1,1); plotmf(fis,'input',1); title('e_n');
subplot(2,1,2); plotmf(fis,'input',2); title('de_n');

%% ---------------- Quick closed-loop test (MATLAB) ----------------
% Plant: 1/(tau s + 1) with ZOH, unity gain (PLC mapping cancels the pump curve)
% One-sample PLC delay on CV. CV clamped [0,100].
T   = 60;  N = round(T/P.Ts);  t = (0:N-1)*P.Ts;
r   = 40*ones(1,N);  r(t >= 10) = 50;          % 10 % step at t = 10 s

y_f = 40*ones(1,N); u_f = 40*ones(1,N); uo_f = 40;   % FLC
y_p = 40*ones(1,N); u_p = 40*ones(1,N); uo_p = 40;   % Dahlin PI + rate clamp
for k = 2:N
    % plant update (uses CV from previous sample -> PLC delay)
    y_f(k) = P.a*y_f(k-1) + (1-P.a)*P.K*uo_f;
    y_p(k) = P.a*y_p(k-1) + (1-P.a)*P.K*uo_p;

    % FLC
    e   = r(k) - y_f(k);
    de  = -(y_f(k) - y_f(k-1))/P.Ts;
    en  = max(-1, min(1, P.Se*e));
    den = max(-1, min(1, P.Sde*de));
    du  = P.Kdu*max(-1, min(1, evalfis(fis,[en den],opt)));
    uo_f = u_f(k-1);
    u_f(k) = max(0, min(100, u_f(k-1) + du));

    % Dahlin PI (P on measurement), same actuator clamp
    e   = r(k) - y_p(k);
    du  = -P.Kc*(y_p(k) - y_p(k-1)) + P.Kc*(1-P.a)*e;
    du  = max(-P.Kdu, min(P.Kdu, du));
    uo_p = u_p(k-1);
    u_p(k) = max(0, min(100, u_p(k-1) + du));
end

figure('Name','Closed loop: 40 -> 50 %');
subplot(2,1,1);
plot(t, r, 'k--', t, y_f, 'LineWidth', 1.3); hold on; plot(t, y_p, 'LineWidth', 1.3);
legend('setpoint','FLC (ANFIS)','Dahlin PI + rate clamp','Location','southeast');
ylabel('flow [%]'); grid on;
subplot(2,1,2);
plot(t, u_f, t, u_p, 'LineWidth', 1.3);
legend('CV FLC','CV PI','Location','southeast'); ylabel('CV [%]'); xlabel('t [s]'); grid on;

ts = @(y) t(find(abs(y - 50) > 0.02*10, 1, 'last') + 1) - 10;   % 2 % settling after step
fprintf('--- Step 40->50 ---\n');
fprintf('FLC: overshoot %.2f %%   settling %.1f s\n', max(y_f) - 50, ts(y_f));
fprintf('PI : overshoot %.2f %%   settling %.1f s\n', max(y_p) - 50, ts(y_p));
fprintf('Ideal (rate-limited) ~ %.1f s\n', 10/P.rate + P.tau);

%% ---------------- Save for Simulink ----------------
writeFIS(fis, 'FLC_F1_anfis');
FLC_F1 = fis;
P1 = P;
save('FLC_F1_anfis_params.mat', 'FLC_F1', 'P1');

%% ---------------- Simulink changes (vs your current model) ----------------
% 1) e path:   gain P1.Se            -> Sat [-1 1] -> FIS in1
% 2) de path:  use y (not e):  y -> (y_k - y_{k-1}) -> gain -P1.Sde/P1.Ts -> Sat [-1 1] -> FIS in2
%              (P on measurement: no kick on setpoint steps)
% 3) FIS block: load FLC_F1 (or FLC_F1_anfis.fis), sample time Ts
% 4) Output:   gain P1.Kdu (= 0.1) -> clamp [-0.1 0.1] -> accumulator -> Sat [0 100]
% 5) Optional: z^-1 on CV before B1_map (PLC delay; also breaks the algebraic loop)