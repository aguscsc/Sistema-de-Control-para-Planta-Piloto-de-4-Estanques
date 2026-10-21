%% build_FLC_F1.m
% Flow loop FIC1 (pump 1) - Phase 2: linear-equivalent fuzzy PI (Dahlin)
%
% Controller (velocity form):
%   du_n = FIS( sat(Se*e), sat(Sde*de) )      normalized, in [-1, 1]
%   du   = Kdu * du_n                          [% per sample]
%   u_k  = sat( u_{k-1} + du , 0, 100 )        accumulator with clamp
%
% With uniform triangles (partition of unity), product AND, diagonal rules,
% the FIS surface is EXACTLY the plane  du_n = 0.5*(e_n + de_n),
% so the FLC reproduces the Dahlin PI:
%   du = Kc*De + Kc*(1-a)*e
%
% Sign convention: e = r - y (positive error -> increase pump).
% de = -(y_k - y_{k-1})/Ts  (P on measurement, no kick on setpoint steps)
% or  (e_k - e_{k-1})/Ts    (classic, kicks on setpoint steps)



%% ---------------- Plant / design parameters ----------------
P.Ts     = 0.1;          % [s] PLC sample time
P.K      = 1.0;          % plant gain seen by controller (mapping cancels pump curve)
P.tau    = 1.3776;       % [s] slowest segment (worst case for design)
P.rate   = 1.0;          % [%/s] pump rate limit
P.lambda = P.tau/2;      % [s] desired closed-loop time constant  <-- design knob

FIS_TYPE = 'sugeno';    % 'mamdani' or 'sugeno' (both give the exact plane)

%% ---------------- Dahlin PI ----------------
P.a  = exp(-P.Ts/P.tau);                 % plant pole (ZOH)
P.q  = exp(-P.Ts/P.lambda);              % desired closed-loop pole
P.Kc = (1 - P.q) / (P.K*(1 - P.a));      % proportional gain
P.Ti = P.Ts / (1 - P.a);                 % integral time (~ tau)

%% ---------------- Scaling factors ----------------
% Output scaling fixed by the actuator: full FIS output = max move per sample
P.Kdu = P.rate * P.Ts;                   % [%/sample] = 0.1

% Matching  Kdu*0.5*(Se*e + Sde*De/Ts)  with  Kc*De + Kc*(1-a)*e
P.Sde = 2*P.Kc*P.Ts / P.Kdu;             % [1/(%/s)]
P.Se  = 2*P.Kc*(1 - P.a) / P.Kdu;        % [1/%]

P.e_sat  = 1/P.Se;                       % [%]   error where e_n saturates
P.de_sat = 1/P.Sde;                      % [%/s] rate where de_n saturates

fprintf('--- Dahlin PI (lambda = %.3f s) ---\n', P.lambda);
fprintf('a = %.4f   q = %.4f   Kc = %.4f   Ti = %.4f s\n', P.a, P.q, P.Kc, P.Ti);
fprintf('--- Scaling factors ---\n');
fprintf('Kdu = %.4f %%/sample   Se = %.4f 1/%%   Sde = %.4f 1/(%%/s)\n', P.Kdu, P.Se, P.Sde);
fprintf('e saturates at +/-%.3f %%,  de saturates at +/-%.3f %%/s\n', P.e_sat, P.de_sat);

%% ---------------- FIS ----------------
switch FIS_TYPE
    case 'mamdani'
        fis = mamfis('Name','FLC_F1', ...
            'AndMethod','prod', 'OrMethod','probor', ...
            'ImplicationMethod','prod', 'AggregationMethod','sum', ...
            'DefuzzificationMethod','centroid');
    case 'sugeno'
        fis = sugfis('Name','FLC_F1', ...
            'AndMethod','prod', 'OrMethod','probor', ...
            'DefuzzificationMethod','wtaver');
end

% Inputs: 3 uniform triangles, partition of unity on [-1, 1]
% (inputs are saturated to [-1,1] before the FIS, so N and P act as shoulders)
fis = addInput(fis, [-1 1], 'Name', 'e');
fis = addInput(fis, [-1 1], 'Name', 'de');
for v = ["e","de"]
    fis = addMF(fis, v, 'trimf', [-2 -1 0], 'Name', 'N');
    fis = addMF(fis, v, 'trimf', [-1  0 1], 'Name', 'Z');
    fis = addMF(fis, v, 'trimf', [ 0  1 2], 'Name', 'P');
end

% Output: 5 terms at -1, -0.5, 0, 0.5, 1
outC = [-1 -0.5 0 0.5 1];
outN = ["NB" "NS" "ZE" "PS" "PB"];
if strcmp(FIS_TYPE,'mamdani')
    % equal-width symmetric triangles + prod implication + sum aggregation
    % -> centroid equals the weighted average of the centers (exact plane)
    fis = addOutput(fis, [-1.5 1.5], 'Name', 'du');
    for k = 1:5
        fis = addMF(fis, 'du', 'trimf', outC(k) + [-0.5 0 0.5], 'Name', outN(k));
    end
else
    fis = addOutput(fis, [-1 1], 'Name', 'du');
    for k = 1:5
        fis = addMF(fis, 'du', 'constant', outC(k), 'Name', outN(k));
    end
end

% Rules: diagonal (MacVicar-Whelan) table
%            de:  N     Z     P
R = [ "NB" "NS" "ZE" ;    % e = N
      "NS" "ZE" "PS" ;    % e = Z
      "ZE" "PS" "PB" ];   % e = P
inT = ["N" "Z" "P"];
rules = strings(0);
for i = 1:3
    for j = 1:3
        rules(end+1) = "e==" + inT(i) + " & de==" + inT(j) + " => du=" + R(i,j); %#ok<SAGROW>
    end
end
fis = addRule(fis, rules);

%% ---------------- Checks ----------------
% 1) Planarity: FIS must equal 0.5*(e+de) everywhere
[E, DE] = meshgrid(linspace(-1,1,41));
opt = evalfisOptions('OutOfRangeInputValueMessage','none', ...
                     'NoRuleFiredMessage','none', ...
                     'EmptyOutputFuzzySetMessage','none');
DU = reshape(evalfis(fis, [E(:) DE(:)], opt), size(E));
err = max(abs(DU(:) - 0.5*(E(:) + DE(:))));
fprintf('--- Check ---\nmax |FIS - plane| = %.2e  (should be ~0)\n', err);

% 2) Surface
figure('Name','FLC_F1 surface');
surf(E, DE, DU); xlabel('e_n'); ylabel('de_n'); zlabel('du_n');
title(sprintf('FLC\\_F1 (%s) - linear-equivalent', FIS_TYPE)); grid on;

%% ---------------- Save for Simulink ----------------
writeFIS(fis, 'FLC_F1');      % creates FLC_F1.fis
FLC_F1 = fis;                 % also in workspace for the Fuzzy Logic Controller block
P1 = P;                       % parameters for pump 1 loop
save('FLC_F1_params.mat', 'FLC_F1', 'P1');

%% ---------------- Simulink wiring (reference) ----------------
% Fixed-step discrete solver, Ts = P1.Ts for all blocks.
%
%  r ──►(+)── e ──► [Gain P1.Se] ──► [Sat ±1] ──► FIS in1
%        ▲ (-)
%        y
%  y ──► [ -(1 - z^-1)/Ts ] ──► [Gain P1.Sde] ──► [Sat ±1] ──► FIS in2   (P on measurement)
%
%  FIS out ──► [Gain P1.Kdu] ──► du ──►(+)──► [Sat 0..100] ──┬──► u (CV_F1) ──► mapping g ──► pump
%                                      ▲                     │
%                                      └──── [Unit Delay] ◄──┘   (IC = initial CV)
%
% Saturation INSIDE the accumulator loop = clamp anti-windup.
% |du| <= Kdu = rate*Ts, so the rate limiter should never be active
% (keep it in the plant model anyway as a check).
%
% Baseline to compare (same loop, no FIS):
%   du = P1.Kc*(-(y_k - y_{k-1})) + P1.Kc*(1 - P1.a)*e_k