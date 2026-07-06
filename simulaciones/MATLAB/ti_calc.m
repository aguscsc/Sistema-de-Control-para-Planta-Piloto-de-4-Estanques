A = 0.159;
a3 = 3.5752e-4;
a4 = 4.1820e-4;
g = 9.81;

h_vals = [0.1, 0.2, 0.3, 0.4, 0.5];
fprintf('h0(m)   T3(s)    T4(s)\n');
for h0 = h_vals
    T3 = (A/a3) * sqrt(2*h0/g);
    T4 = (A/a4) * sqrt(2*h0/g);
    fprintf('%.1f    %.1f    %.1f\n', h0, T3, T4);
end
