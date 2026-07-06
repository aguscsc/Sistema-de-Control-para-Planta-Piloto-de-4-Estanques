% Busca, entre las 127 semillas no nulas del LFSR de 7 bits (x^7+x^6+1,
% igual al de horext_prof_test.m: s=sh(7)+sh(6), XOR via umbral 1.5),
% la semilla que maximiza la cantidad de transiciones (cambios 0<->1)
% en las primeras N salidas. Esto evita arrancar el PRBS metido en una
% corrida larga del mismo simbolo.

N = 20;   % ventana inicial a evaluar

best_count   = -1;
best_seed    = [];
best_outputs = [];

for seed_int = 1:127
    sh = double(bitget(seed_int, 7:-1:1));   % semilla de 7 bits (sh1..sh7)

    out = zeros(1, N);
    for k = 1:N
        s = sh(7) + sh(6);
        if s >= 1.5
            s = 0;
        end
        sh = [s, sh(1:6)];   % shift: nuevo bit entra en sh1, sh7 se descarta
        out(k) = s;
    end

    transitions = sum(abs(diff(out)));

    if transitions > best_count
        best_count   = transitions;
        best_seed    = sh; %#ok<NASGU> % (se sobreescribe, guardamos la semilla original abajo)
        best_seed    = double(bitget(seed_int, 7:-1:1));
        best_outputs = out;
    end
end

fprintf('Mejor semilla (sh1..sh7): %s\n', mat2str(best_seed));
fprintf('Transiciones en las primeras %d salidas: %d\n', N, best_count);
fprintf('Salidas: %s\n', mat2str(best_outputs));
