% =========================================================================
% Parámetros de la Planta de 4 Estanques (Sistema Internacional)
% Hardware: Danfoss VLT 2800 + Bomba Shertec Hypro (1 HP)
% =========================================================================
clc
clear
% Constante de Gravedad
g = 9.81; % [m/s^2]

% Áreas transversales de los estanques (Ai)
% Estimación para estanques cilíndricos de laboratorio
A = 0.159; A2 = 0.159; A3 = 0.159; A4 = 0.159; % [m^2]

% Áreas de las válvulas de descarga/fugas (ai)
a1 =  6.1360e-04; a2 = 5.4678e-04; a3 = 3.5752e-04; a4 = 4.1820e-04; % [m^2]
%a1 =  0.00497; a2 = 0.00497; a3 = 0.00497; a4 =  0.00497;
% Posición de las válvulas de 3 vías (Gamma)
gamma1 = 0.455; % gamma1% del flujo de Bomba 1 va al Estanque 1
gamma2 = 0.4195; % gamma2% del flujo de Bomba 2 va al Estanque 2

q_max = 0.0017; % caudal máximo supuesto (m^3 / s)
h_max = 0.6; % estanque de 60cm
k1 = q_max;
k2 = q_max;
%% calculo de limitaciones

% max altura acción individual estanque 1
syms h1 h2 h3 h4
eq1 = gamma1*k1*1 - a1*sqrt(2*g*h1);
h1_max_alone = double(solve(eq1,h1));
eq1_2 = (gamma1*k1*1) + a3*sqrt(2*g*h3)== a1*sqrt(2*g*h_max);
h3_for_h1 = double(solve(eq1_2,h3));
fprintf('El Estanque 1 necesita que h3 = %.3f metros\n', h3_for_h1);

eq3 = (1-gamma2)*k2*1 - a3*sqrt(2*g*h3);
h3_max = double(solve(eq3,h3));
fprintf('El Estanque 3 alcanza %.3f metros\n', h3_max);

eq2 = gamma2*k2*1 - a2*sqrt(2*g*h2);
h2_max_alone = double(solve(eq2,h2));
eq2_2 = (gamma2*k2*1) + a4*sqrt(2*g*h4)== a2*sqrt(2*g*h_max);
h4_for_h2 = double(solve(eq2_2,h4));
fprintf('El Estanque 2 necesita que h4 = %.3f metros\n', h4_for_h2);


eq4 = (1-gamma1)*k1*1 - a4*sqrt(2*g*h4);
h4_max = double(solve(eq4,h4));
fprintf('El Estanque 4 alcanza %.3f metros\n', h4_max);

%% max h
eqmax1 = gamma1*k1*1 - a1*sqrt(2*g*h1) + a3*sqrt(2*g*h3_max);
h1_max = double(solve(eqmax1, h1));
fprintf('El estanque 1 alcaza %.3f metros, %.3f\n', h1_max, h1_max*(100/h_max))

eqmax2 = gamma2*k2*1 - a2*sqrt(2*g*h2) + a4*sqrt(2*g*h4_max);
h2_max = double(solve(eqmax2, h2));
fprintf('El estanque 2 alcaza %.3f metros, %.3f\n', h2_max, h2_max*(100/h_max))

fprintf('El estanque 3 alcaza %.3f metros, %.3f \n', h3_max, h3_max*(100/h_max))
fprintf('El estanque 4 alcaza %.3f metros, %.3f \n', h4_max, h4_max*(100/h_max))

%% prueba 60%

n4 = ((1-gamma1)*k1*0.6)/(a4*sqrt(2*g));
n4r = n4^2
n3 = ((1-gamma2)*k2*0.6)/(a3*sqrt(2*g));
n3r = n3^2
n2 = (a4*sqrt(2*g)*n4 + (gamma2)*k2*0.6)/(a2*sqrt(2*g));
n2r = n2^2
n1 = (a3*sqrt(2*g)*n3 + (gamma1)*k1*0.6)/(a1*sqrt(2*g));
n1r = n1^2