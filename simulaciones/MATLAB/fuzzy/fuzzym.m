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