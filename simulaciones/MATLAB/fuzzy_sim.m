%% quadtank_sim.m  --  Johansson 4-tank plant, pure MATLAB (no toolboxes)
% Proyecto 1 MFLC.  Purpose: debug harness for the fuzzy controllers.
%
%   Tanks 1,2 = lower ("master") tanks, tanks 3,4 = upper tanks.
%   Pump 1 -> tank 1 (gamma1) and tank 4 (1-gamma1)
%   Pump 2 -> tank 2 (gamma2) and tank 3 (1-gamma2)
%   Tank 3 drains into tank 1, tank 4 drains into tank 2.
%   Model: main.pdf eqs. (1)-(4), Torricelli outflow, u in [0,1], q_pump = k*u.
%
% Needs R2016b+ (local functions at the end of a script file).
% Controller interface:  [u, cs] = ctrl(kk, t, y, r, cs, P)
%   y  = 4x1 measured levels [h1;h2;h3;h4]   (h3,h4 available for the FF)
%   r  = 2x1 level references for tanks 1,2
%   u  = 2x1 pump commands in [0,1] (clipped again by the harness)
%   cs = controller state struct (persisted by the harness)

clear; clc; close all;

%% 1) Parameters  (main.pdf Table 1.1; gamma2 = 0.4195 as in sec. 7.4)
P.g     = 9.81;                                   % [m/s^2]
P.A     = 0.159*ones(1,4);                        % tank areas [m^2]
P.a     = [6.1360e-04 5.4678e-04 3.5752e-04 4.1820e-04];   % orifice areas [m^2]
P.gamma = [0.455 0.4195];                         % 3-way valve splits
P.k     = [0.0017 0.0017];                        % pump gain at u=1 [m^3/s]
P.hmax  = 0.6;                                    % tank height [m] (overflow clip)
P.Ts    = 0.1;                                    % sample time [s]
P.nsub  = 10;                                     % RK4 sub-steps per sample

%% 2) Operating point, feasibility, linearization
h_op = [0.20 0.20];                               % lower-tank levels at the OP [m]
OP   = op_point(P, h_op);
fprintf('Operating point h1*=%.3f h2*=%.3f -> u*=[%.3f %.3f], h3*=%.3f h4*=%.3f  (%s)\n', ...
    h_op(1), h_op(2), OP.u(1), OP.u(2), OP.h(3), OP.h(4), ...
    ternary(OP.feasible, 'feasible', 'INFEASIBLE: pump or tank limit violated'));

T  = (P.A ./ P.a) .* sqrt(2*OP.h.' / P.g);        % linearized time constants [s]
fprintf('Time constants T1..T4 = %.1f  %.1f  %.1f  %.1f s\n', T);

Alin = [-1/T(1) 0 P.A(3)/(P.A(1)*T(3)) 0;
        0 -1/T(2) 0 P.A(4)/(P.A(2)*T(4));
        0 0 -1/T(3) 0;
        0 0 0 -1/T(4)];
Blin = [P.gamma(1)*P.k(1)/P.A(1) 0;
        0 P.gamma(2)*P.k(2)/P.A(2);
        0 (1-P.gamma(2))*P.k(2)/P.A(3);
        (1-P.gamma(1))*P.k(1)/P.A(4) 0];
Clin = [1 0 0 0; 0 1 0 0];
G0   = -Clin * (Alin \ Blin);                     % DC gain [m per unit u]
RGA  = G0 .* inv(G0).';
fprintf('DC gain G(0) =\n'); disp(G0);
fprintf('RGA =\n');         disp(RGA);
fprintf('gamma1+gamma2 = %.4f -> %s\n', sum(P.gamma), ...
    ternary(sum(P.gamma) < 1, 'RHP transmission zero (non-minimum phase)', 'minimum phase'));
z = roots([P.gamma(1)*P.gamma(2)*T(3)*T(4), P.gamma(1)*P.gamma(2)*(T(3)+T(4)), sum(P.gamma)-1]);
fprintf('Transmission zeros (eq. 32, equal areas): %s\n', mat2str(z.', 4));

%% 3) Feasible setpoint map (pump limits 0<=u<=1, tank limits)
do_map = true;
if do_map
    hh = 0.10:0.005:0.35;
    Fm = false(numel(hh), numel(hh));              % rows: h2, cols: h1
    for i = 1:numel(hh)
        for j = 1:numel(hh)
            o = op_point(P, [hh(i) hh(j)]);
            Fm(j,i) = o.feasible;
        end
    end
    figure('Name','Feasible setpoints'); imagesc(hh, hh, double(Fm)); axis xy; hold on;
    plot(h_op(1), h_op(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    colormap([0.85 0.85 0.85; 0.55 0.8 0.55]);
    xlabel('h_1^* [m]'); ylabel('h_2^* [m]'); title('Reachable (h_1^*, h_2^*) - green');
end

%% 4) Simulation options
opt.noise_std = 0.0;                              % level sensor noise [m] (try 0.001)
opt.pump_tau  = 0.0;                              % pump/flow-loop lag [s]; 0 = ideal
opt.u_init    = OP.u;

% Reference profile (must stay inside the feasible region from section 3)
ref_fun = @(t) h_op(:) + [0.02*(t >= 50); 0.02*(t >= 500)];
Tend    = 1000;                                   % [s]
h0      = OP.h;                                   % start at equilibrium

%% 5) Scenario A: hold u* -> levels must not move (plant/equilibrium consistency)
csA = struct('u0', OP.u);
RA  = run_sim(P, @ctrl_hold, csA, h0, @(t) h_op(:), 600, opt);
driftA = max(max(abs(RA.h - OP.h.')));
fprintf('\n[A] hold u*: max level drift = %.2e m  (should be ~0)\n', driftA);

%% 6) Scenario B: open-loop step of +0.05 on each pump command (one at a time)
csB1 = struct('u0', OP.u + [0.05; 0]);
csB2 = struct('u0', OP.u + [0; 0.05]);
RB1  = run_sim(P, @ctrl_hold, csB1, h0, @(t) h_op(:), 800, opt);
RB2  = run_sim(P, @ctrl_hold, csB2, h0, @(t) h_op(:), 800, opt);
figure('Name','Open-loop steps');
subplot(1,2,1); plot(RB1.t, RB1.h - OP.h.', 'LineWidth', 1.2); grid on;
title('+0.05 on u_1'); xlabel('t [s]'); ylabel('\Delta h [m]'); legend('h1','h2','h3','h4','Location','best');
subplot(1,2,2); plot(RB2.t, RB2.h - OP.h.', 'LineWidth', 1.2); grid on;
title('+0.05 on u_2'); xlabel('t [s]');

%% 7) Scenario C: closed loop baseline (PI + static FF on the reference)
% pair = [1 2]: level i -> pump i (diagonal, RGA(1,1)<0, gain ceiling)
% pair = [2 1]: h1 -> pump 2, h2 -> pump 1 (off-diagonal)
csC.M    = OP.M;
csC.pair = [1 2];
csC.Kp   = 0.5;                                   % [unit u / m]
csC.Ki   = csC.Kp / 100;                          % Ti = 100 s
csC.I    = zeros(2,1);
RC = run_sim(P, @ctrl_pi, csC, h0, ref_fun, Tend, opt);
report_run('PI baseline', RC, P);
plot_run(RC, P, 'Closed loop: PI baseline (no dynamic FF from upper tanks)');

%% 8) Where the fuzzy blocks plug in
% Copy ctrl_pi to ctrl_fuzzy.  Keep u_ff = M\Fref (static part), replace the PI by
% F11/F22 (e,de) and add the FF block F12/F21 with inputs (h_upper, dh_upper):
%   dh3 = filtered derivative of y(3);  dh4 = filtered derivative of y(4).
% Compare with RC using report_run / plot_run.

%% ===================== local functions =====================

function OP = op_point(P, h12)
% Static equilibrium for given lower-tank levels (main.pdf sec. 8, eqs. 35-42)
g = P.g; a = P.a; gm = P.gamma; k = P.k;
M = [ gm(1)*k(1),      (1-gm(2))*k(2);
     (1-gm(1))*k(1),    gm(2)*k(2)   ];
F = [a(1)*sqrt(2*g*h12(1)); a(2)*sqrt(2*g*h12(2))];
u = M \ F;
h3 = ((1-gm(2))*k(2)*u(2) / (a(3)*sqrt(2*g)))^2;
h4 = ((1-gm(1))*k(1)*u(1) / (a(4)*sqrt(2*g)))^2;
OP.M = M; OP.F = F; OP.u = u;
OP.h = [h12(1); h12(2); h3; h4];
OP.feasible = all(u >= 0) && all(u <= 1) && all(OP.h <= P.hmax);
end

function dh = plant_rhs(h, u, P)
% Nonlinear plant, main.pdf eqs. (1)-(4).  h: 4x1, u: 2x1
h  = max(h, 0);
q  = P.a(:) .* sqrt(2*P.g*h);                     % orifice outflows [m^3/s]
qp = P.k(:) .* u;                                 % pump flows [m^3/s]
dh = [ -q(1) + q(3) + P.gamma(1)*qp(1);
       -q(2) + q(4) + P.gamma(2)*qp(2);
       -q(3) + (1-P.gamma(2))*qp(2);
       -q(4) + (1-P.gamma(1))*qp(1) ] ./ P.A(:);
end

function h = plant_step(h, u, P)
% One sample of the plant with zero-order-hold input, RK4 inside
dt = P.Ts / P.nsub;
for i = 1:P.nsub
    k1 = plant_rhs(h,           u, P);
    k2 = plant_rhs(h + dt/2*k1, u, P);
    k3 = plant_rhs(h + dt/2*k2, u, P);
    k4 = plant_rhs(h + dt*k3,   u, P);
    h  = h + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    h  = min(max(h, 0), P.hmax);                  % empty / overflow limits
end
end

function R = run_sim(P, ctrl, cs, h0, ref_fun, Tend, opt)
% Fixed-step harness: measure -> controller -> saturate -> (pump lag) -> plant
N   = round(Tend / P.Ts);
R.t = (0:N-1).' * P.Ts;
R.h = zeros(N,4); R.y = zeros(N,4); R.u = zeros(N,2); R.r = zeros(N,2);
h   = h0(:);
up  = opt.u_init(:);
for kk = 1:N
    t = R.t(kk);
    r = ref_fun(t);
    y = h + opt.noise_std * randn(4,1);
    [uc, cs] = ctrl(kk, t, y, r(:), cs, P);
    uc = min(max(uc(:), 0), 1);
    if opt.pump_tau > 0
        up = up + P.Ts/(opt.pump_tau + P.Ts) * (uc - up);
    else
        up = uc;
    end
    R.h(kk,:) = h.'; R.y(kk,:) = y.'; R.u(kk,:) = up.'; R.r(kk,:) = r(:).';
    h = plant_step(h, up, P);
end
end

function [u, cs] = ctrl_hold(~, ~, ~, ~, cs, ~)
u = cs.u0;                                        % constant pump commands
end

function [u, cs] = ctrl_pi(~, ~, y, r, cs, P)
% Baseline: static FF from the reference + PI per level loop, conditional anti-windup
e    = r - y(1:2);
cs.I = cs.I + cs.Ki .* e * P.Ts;
a12  = P.a(1:2).';
Fref = a12 .* sqrt(2*P.g*r);                      % outflow needed at the reference
u    = cs.M \ Fref;                               % static FF
du   = cs.Kp .* e + cs.I;
for i = 1:2
    u(cs.pair(i)) = u(cs.pair(i)) + du(i);
end
us = min(max(u, 0), 1);
if any(abs(u - us) > 0)
    cs.I = cs.I - cs.Ki .* e * P.Ts;              % freeze integrator while saturated
end
u = us;
end

function report_run(name, R, P)
e  = R.h(:,1:2) - R.r;
fprintf('\n[%s]\n', name);
fprintf('  IAE h1 = %.3f  h2 = %.3f  [m*s]\n', sum(abs(e(:,1)))*P.Ts, sum(abs(e(:,2)))*P.Ts);
fprintf('  final error h1 = %+.4f  h2 = %+.4f [m]\n', e(end,1), e(end,2));
fprintf('  pump command range u1 = [%.2f %.2f]  u2 = [%.2f %.2f]\n', ...
    min(R.u(:,1)), max(R.u(:,1)), min(R.u(:,2)), max(R.u(:,2)));
fprintf('  time at saturation: u1 %.1f%%  u2 %.1f%%\n', ...
    100*mean(R.u(:,1) >= 0.999 | R.u(:,1) <= 0.001), 100*mean(R.u(:,2) >= 0.999 | R.u(:,2) <= 0.001));
end

function plot_run(R, P, ttl)
hs = movmean(R.h(:,3:4), 30, 1);                  % smooth, then differentiate in time
dh = [diff(hs, 1, 1); zeros(1,2)] / P.Ts;         % rate of the upper tanks [m/s]
figure('Name', ttl);
subplot(4,1,1); plot(R.t, R.h(:,1:2), R.t, R.r, '--'); grid on;
ylabel('h_{1,2} [m]'); title(ttl); legend('h1','h2','r1','r2','Location','best');
subplot(4,1,2); plot(R.t, R.h(:,3:4)); grid on; ylabel('h_{3,4} [m]'); legend('h3','h4','Location','best');
subplot(4,1,3); plot(R.t, dh*1e3); grid on; ylabel('dh_{3,4}/dt [mm/s]'); legend('dh3','dh4','Location','best');
subplot(4,1,4); plot(R.t, R.u); grid on; ylabel('u'); xlabel('t [s]'); legend('u1','u2','Location','best'); ylim([-0.05 1.05]);
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end