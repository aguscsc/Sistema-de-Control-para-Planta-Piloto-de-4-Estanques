% Calculo de parametros PRBS para identificacion del lazo de nivel (L3)
% Misma fisica que ti_calc.m (Torricelli, T_i = (A/a)*sqrt(2*h0/g)),
% aplicada ahora al diseno del test: Ts, delta y duracion deben escalar
% con la constante de tiempo real del tanque (~80-100x mas lenta que F1).

A  = 0.159;
a3 = 3.5752e-4;
g  = 9.81;

% --- Punto de operacion de referencia ---
% Ajustar h0_ref al nivel real en torno al cual se hara el PRBS
h0_ref = 0.3;   % m

tau3 = (A/a3) * sqrt(2*h0_ref/g);
fprintf('tau3 en h0=%.2f m   = %.1f s\n', h0_ref, tau3);

% --- Ts: ~5 muestras por constante de tiempo (regla estandar de ID) ---
Ts = round(tau3/5/5)*5;        % redondeado a multiplo de 5 s
p3 = exp(-Ts/tau3);
fprintf('Ts recomendado       = %d s   (polo discreto p = %.3f)\n', Ts, p3);

% --- Delta: 3.5x tau, para evitar colinealidad rampa/delta (ver caso F1) ---
delta = round(3.5*tau3/Ts)*Ts;
fprintf('Delta recomendado    = %d s   (%d muestras a Ts=%ds)\n', delta, delta/Ts, Ts);

% --- Duracion: 10-12 periodos delta, suficientes conmutaciones para RLS ---
dur_min = 10*delta;
dur_max = 12*delta;
fprintf('Duracion recomendada = %d - %d s  (~%.1f - %.1f min)\n', ...
    dur_min, dur_max, dur_min/60, dur_max/60);

% --- Sensibilidad de tau a la amplitud (para acotar el swing de nivel) ---
dh = 0.08;  % +/- 8 cm de excursion alrededor de h0_ref
tau3_lo = (A/a3)*sqrt(2*(h0_ref-dh)/g);
tau3_hi = (A/a3)*sqrt(2*(h0_ref+dh)/g);
fprintf('tau3 en h0=%.2f m   = %.1f s  (limite inferior del swing)\n', h0_ref-dh, tau3_lo);
fprintf('tau3 en h0=%.2f m   = %.1f s  (limite superior del swing)\n', h0_ref+dh, tau3_hi);

% --- Caudal de entrada en estado estacionario para sostener h0_ref ---
qss = a3*sqrt(2*g*h0_ref);     % m^3/s
fprintf('Caudal estacionario en h0=%.2f m = %.4f L/s (%.2f L/min)\n', ...
    h0_ref, qss*1000, qss*1000*60);

% --- a0 teorico = p + p^2 + ... + p^T, para T=1,2,3 (con el Ts elegido) ---
fprintf('\n--- a0 teorico (p=%.4f a Ts=%d s) ---\n', p3, Ts);
suma = 0;
for T = 1:3
    suma = suma + p3^T;
    fprintf('T=%d : a0 = %.4f\n', T, suma);
end
