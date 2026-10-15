%% deadbeat_flc_pump2.m
% Builds a 3x3 incremental fuzzy controller (inputs e, de; output du) that
% behaves like a deadbeat / Dahlin controller for pump 2, using the
% Takagi-Sugeno one-step model fis_tuned_flujo_2f.fis.
%
% Steps:
%   1. Identify K and tau (in samples) from step tests on the TS model,
%      keeping only operating points where the model is consistent.
%   2. Find u_fold: the pump % above which the model stops being
%      monotonic in u (more pump -> less flow). Not physical.
%   3. Deadbeat/Dahlin design for q(k+1) = a q(k) + b u(k):
%         du(k) = (1-beta)/b * ( e(k) - a e(k-1) )
%      -> incremental PI, mapped onto GE, GdE, GdU.
%   4. Build the FLC (prod/probor, 3 triangular MFs, singletons) and save it.
%   5. Simulate it on the linear surrogate and on the TS model.
%
% Assumes the model is: in1 = q(k) [%], in2 = u(k) [%], out = q(k+1) [%].
% The FLC must run at the SAME sample time the model was trained with.
%
% Implementation in Simulink:
%   E  = sat(GE*e, -1, 1),  dE = sat(GdE*(e(k)-e(k-1)), -1, 1)
%   ui(k) = sat( ui(k-1) + GdU*FLC(E,dE) - kb*(ui(k-1) - u(k-1)), 0, 100 )
%   u(k)  = u(k-1) + sat( ui(k) - u(k-1), -rateMax*Ts, rateMax*Ts )
%   (ui = controller's internal command, u = what is sent to the pump)

clear; clc; close all;

%% ---------------- Parameters ----------------
modelFile = 'fis_tuned_flujo_2f.fis';
outFile   = 'FLC_deadbeat_pump2.fis';

lambda  = 2;      % desired closed-loop time constant [samples]; 0 = pure deadbeat
stepLin = 10;      % largest setpoint step [%] the FLC should treat linearly
sp0 = 40;       % test: initial setpoint [%]  (keep inside the clean region)
sp1 = 50;         % test: final setpoint [%]
N   = 300;        % simulation length [samples]

% Actuator limits
Ts      = 1;      % sample time [s] -- MUST be the one the model was trained with
rateMax = 1;      % pump slew-rate limit [%/s]
kb      = 0.05;   % back-calculation gain [0..1]: 0 = no anti-windup (overshoot),
                  % 1 = clamp du directly (no overshoot, slow). Tune 0.03-0.15.

uOP     = (35:5:85)';  % operating points for identification [% pump]
duId    = 3;           % identification step size [% pump]
Nid     = 3000;        % samples per identification run
qCheck  = [20 60];     % flow range [%] where monotonicity in u is checked

%% ---------------- Load model ----------------
opt = evalfisOptions('OutOfRangeInputValueMessage','none', ...
                     'NoRuleFiredMessage','none', ...
                     'EmptyOutputFuzzySetMessage','none');
model = readfis(modelFile);
f = @(q,u) evalfis(model, [q(:) u(:)], opt);   % vectorized: q(k+1) = f(q(k),u(k))

%% ---------------- 1. Identify K and tau ----------------
nOP = numel(uOP);
q = 0.8*uOP;
for k = 1:Nid, q = f(q, uOP); end            % settle at each operating point
q0 = q;

s    = [duId -duId];
Kid  = nan(nOP,2);  tauId = nan(nOP,2);  mono = false(nOP,2);
for j = 1:2
    q  = q0;  tr = zeros(Nid, nOP);
    for k = 1:Nid
        q = f(q, uOP + s(j));
        tr(k,:) = q';
    end
    dq = tr(end,:)' - q0;
    Kid(:,j) = dq / s(j);
    for i = 1:nOP
        idx = find(abs(tr(:,i) - q0(i)) >= 0.632*abs(dq(i)), 1);
        if ~isempty(idx), tauId(i,j) = idx; end
        mono(i,j) = all(diff(tr(:,i))*sign(dq(i)) >= -1e-9);
    end
end

% keep points where +step and -step agree and the response is sensible
ok = all(Kid > 0, 2) & all(mono, 2) & ...
     abs(Kid(:,1)-Kid(:,2))     ./ max(Kid,[],2)   < 0.15 & ...
     abs(tauId(:,1)-tauId(:,2)) ./ max(tauId,[],2) < 0.30;

disp('Identification on the TS model (tau in samples):');
disp(table(uOP, q0, Kid(:,1), Kid(:,2), tauId(:,1), tauId(:,2), ok, ...
     'VariableNames', {'u0','q0','K_up','K_down','tau_up','tau_down','used'}));

Km   = median(Kid(ok,:),   'all');
tauM = median(tauId(ok,:), 'all');
a = exp(-1/tauM);
b = Km*(1 - a);
fprintf('Plant: K = %.3f, tau = %.1f samples  ->  a = %.4f, b = %.4f\n', Km, tauM, a, b);

%% ---------------- 2. Fold check ----------------
uGrid = linspace(0, 100, 1001)';
uFold = 100;
for qq = qCheck(1):1:qCheck(2)
    v = f(qq*ones(size(uGrid)), uGrid);
    i = find(diff(v) <= 0, 1);
    if ~isempty(i), uFold = min(uFold, uGrid(i)); end
end
fprintf('Model is monotonic in u only below u_fold = %.1f %%\n', uFold);

%% ---------------- 3. Deadbeat / Dahlin gains ----------------
beta = 0;
if lambda > 0, beta = exp(-1/lambda); end
Kp   = (1-beta)*a/b;          % proportional gain
KiTs = (1-beta)*(1-a)/b;      % integral gain per sample

% Near the origin the FLC is F = E + dE, so:
%   Kp = GdU*GdE,  KiTs = GdU*GE   ->  GdE/GE = a/(1-a)  (fixed by the plant)
% eMax chosen so that dE reaches 1 exactly when de = stepLin.
eMax = stepLin*a/(1-a);
GE   = 1/eMax;
GdE  = GE*a/(1-a);
GdU  = KiTs/GE;

fprintf('\nlambda = %g samples  ->  Kp = %.2f, Ki*Ts = %.3f\n', lambda, Kp, KiTs);
fprintf('FLC gains: GE = %.5f   GdE = %.4f   GdU = %.2f\n', GE, GdE, GdU);
fprintf('(e saturates at +-%.0f %%, de saturates at +-%.1f %% per sample)\n', eMax, 1/GdE);

%% ---------------- 4. Build the FLC ----------------
flc = sugfis('Name', 'FLC_deadbeat_pump2', ...
    'AndMethod', 'prod', 'OrMethod', 'probor', ...
    'ImplicationMethod', 'prod', 'AggregationMethod', 'sum', ...
    'DefuzzificationMethod', 'wtaver');

inNames = {'e', 'de'};
for c = 1:2
    nm  = inNames{c};
    flc = addInput(flc, [-1 1], 'Name', nm);
    flc = addMF(flc, nm, 'trapmf', [-100 -100 -1 0], 'Name', 'N');
    flc = addMF(flc, nm, 'trimf',  [-1 0 1],         'Name', 'Z');
    flc = addMF(flc, nm, 'trapmf', [0 1 100 100],    'Name', 'P');
end
flc = addOutput(flc, [-1 1], 'Name', 'du');
flc = addMF(flc, 'du', 'constant', -1, 'Name', 'N');
flc = addMF(flc, 'du', 'constant',  0, 'Name', 'Z');
flc = addMF(flc, 'du', 'constant',  1, 'Name', 'P');

% Rule table R(e,de): 1=N, 2=Z, 3=P
R = [1 1 2;     % e = N
     1 2 3;     % e = Z
     2 3 3];    % e = P
rules = zeros(9,5);  r = 0;
for i = 1:3
    for j = 1:3
        r = r + 1;
        rules(r,:) = [i j R(i,j) 1 1];   % [e de du weight AND]
    end
end
flc = addRule(flc, rules);
writeFIS(flc, outFile);
fprintf('\nSaved %s\n', outFile);

% sanity check: in the linear quadrant F must equal E + dE
fprintf('Check: F(0.3,-0.2) = %.4f (should be 0.1000)\n', evalfis(flc, [0.3 -0.2], opt));

%% ---------------- 5. Closed-loop simulation ----------------
plants = {@(q,u) a*q + b*u,  f};
pNames = {'Linear surrogate (K, \tau)', sprintf('TS model, u \\leq %.1f%%', uFold)};
uMax   = [100, uFold];

figure('Name', 'Deadbeat-like FLC, pump 2', 'Color', 'w');
for p = 1:2
    % start at equilibrium for sp0
    if p == 1
        u = sp0/Km;  q = sp0;
    else
        ug = linspace(0, uFold, 2001)';
        [~, i] = min(abs(f(sp0*ones(size(ug)), ug) - sp0));
        u = ug(i);  q = sp0;
        for k = 1:3000, q = f(q, u); end
    end

    dMax  = rateMax*Ts;        % max pump change per sample [%]
    ui    = u;                 % controller's internal (desired) pump command
    ePrev = sp0 - q;
    Q = zeros(N,1);  U = zeros(N,1);  UI = zeros(N,1);  SP = zeros(N,1);
    for k = 1:N
        sp = sp0 + (sp1 - sp0)*(k > 10);
        e  = sp - q;
        de = e - ePrev;  ePrev = e;
        E  = max(min(GE*e,  1), -1);
        dE = max(min(GdE*de, 1), -1);
        F  = evalfis(flc, [E dE], opt);

        % accumulator with back-calculation anti-windup
        ui = ui + GdU*F - kb*(ui - u);
        ui = min(max(ui, 0), uMax(p));
        % actuator: slew-rate limit, then range limit
        u  = u + min(max(ui - u, -dMax), dMax);
        u  = min(max(u, 0), uMax(p));

        q  = plants{p}(q, u);
        Q(k) = q;  U(k) = u;  UI(k) = ui;  SP(k) = sp;
    end

    kOut = find(abs(Q(11:end) - sp1) >= 0.2, 1, 'last');
    if isempty(kOut), kOut = 0; end
    fprintf('%-28s settles (+-0.2%%) in %3d samples (%.1f s), overshoot %.2f %%\n', ...
            strrep(pNames{p}, '\', ''), kOut, kOut*Ts, max(Q(11:end)) - sp1);

    t = (0:N-1)'*Ts;
    subplot(2,2,p);
    stairs(t, SP, 'k--'); hold on; stairs(t, Q, 'LineWidth', 1.4); grid on;
    ylabel('flow [%]'); title(pNames{p}); legend('SP', 'q', 'Location', 'southeast');
    subplot(2,2,p+2);
    stairs(t, UI, ':', 'LineWidth', 1.2); hold on; stairs(t, U, 'LineWidth', 1.4); grid on;
    ylabel('pump [%]'); xlabel('time [s]'); legend('u_{int}', 'u applied', 'Location', 'southeast');
end