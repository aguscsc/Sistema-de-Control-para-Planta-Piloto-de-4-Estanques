clc
clear
% Constante de Gravedad
g = 9.81; % [m/s^2]

% Posición de las válvulas de 3 vías (Gamma)
gamma1 = 0.455; % gamma1% del flujo de Bomba 1 va al Estanque 1
gamma2 = 0.4195; % gamma2% del flujo de Bomba 2 va al Estanque 2

q_max = 0.0017; % caudal máximo supuesto (m^3 / s)
h_max = 0.6; % estanque de 60cm
k1 = q_max;
k2 = q_max;
%% Primera medición
FIT01 = 60*q_max/100;
FIT02 = 60*q_max/100;
h1 = 11.41*h_max/100;
h2 = 5.6*h_max/100;
h3 = 17.6*h_max/100;
h4 = 21.50*h_max/100;
% calculo de areas de descarga
a41 = ((1-gamma1)*FIT01)/(sqrt(2*g*h4))
a31 = ((1-gamma2)*FIT02)/(sqrt(2*g*h3))
a21 = (gamma2*FIT02)/(sqrt(2*g*h2))
a11 = (gamma1*FIT01)/(sqrt(2*g*h1))

%% segunda medición
LIT01 = 60*q_max/100;
LIT02 = 60*q_max/100;
h1 = 25.17*h_max/100;
h2 = 27.5*h_max/100;
h3 = 23.3*h_max/100;
h4 = 15.01*h_max/100;
% calculo de areas de descarga
a42 = ((1-gamma1)*LIT01)/(sqrt(2*g*h4))
a32 = ((1-gamma2)*LIT02)/(sqrt(2*g*h3))
a22 = (a42*sqrt(2*g*h4) + gamma2*LIT02)/(sqrt(2*g*h2))
a12 = (a32*sqrt(2*g*h3) + gamma1*LIT01)/(sqrt(2*g*h1))

%% avg
a1 = (a11 + a12)/2
a2 = mean([a21 a22],2)
a3 = mean([a31 a32],2)
a4 = mean([a41 a42],2)