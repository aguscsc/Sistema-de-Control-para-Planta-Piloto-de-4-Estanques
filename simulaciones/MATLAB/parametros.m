% =========================================================================
% Parámetros de la Planta de 4 Estanques (Sistema Internacional)
% Hardware: Danfoss VLT 2800 + Bomba Shertec Hypro (1 HP)
% =========================================================================

% Constante de Gravedad
g = 9.81; % [m/s^2]

% Áreas transversales de los estanques (Ai)
% Estimación para estanques cilíndricos de laboratorio
A1 = 0.172; A2 = 0.172; A3 = 0.172; A4 = 0.172; % [m^2]

% Áreas de las válvulas de descarga/fugas (ai)
a1 =  6.3347e-4; a2 = 6.3347e-4; a3 = 6.3347e-4; a4 =  6.3347e-4; % [m^2]

% 4. Posición de las válvulas de 3 vías (Gamma)
gamma1 = 0.455; % gamma1% del flujo de Bomba 1 va al Estanque 1
gamma2 = 0.455; % gamma2% del flujo de Bomba 2 va al Estanque 2

q_max = 0.003; % caudal máximo supuesto (m^3 / s)
h_max = 0.2; % estanque de 20cm, falta correción
k1 = 0.00006;
k2 = 0.00006;
max = k1*50;