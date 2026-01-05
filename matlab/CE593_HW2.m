%Parameters
h = 0.4; g = 9.81;
%f = linspace(0.1, 10, 50);
%df = f(2) - f(1);

edges = linspace(0.1, 10, 51);
widths = diff(edges);
df = widths(1);
f = edges(1:end-1) + widths .* rand(1, 50);
N = numel(f);

T = 1 ./ f ;
L = arrayfun(@(t) dispersion_equation(t, h), T);
k = 2 * pi ./ L;

%Equations
n = 0.5 .* (1 + (2 * k * h) ./ (sinh(2 * k * h)));
phi = 0.5 .* tanh(k * h).^2 ./ n;
a = sqrt(2 * df .* (0.205 * 0.15^2 * 1.6^-4 .* f.^-5 .* exp(-0.75 * ...
    (1.6 .* f).^-4)) .* phi);
epsilon = 2 * pi * rand(1,N);

%Time ranges
t1 = 0 : 0.05 : 10;
t2 = 0 : 0.05 : 300;

%Part 1
eta1 = a.' .* cos(2 * pi * f.' .* t1 + epsilon.');

figure(1); hold on; grid on;
for idx = [1 10 30 50]
    plot(t1, eta1(idx,:), 'DisplayName', sprintf('\\eta_{%d}(t)', idx), 'LineWidth', 1.1);
end

legend
title('Individual surface elevations (0–10 s)');
xlabel('Time (s)')
ylabel('\eta(t) (m)')

%Part 2
eta2 = sum(a.' .* cos(2 * pi * f.' .* t2 + epsilon.'), 1);

figure(2); grid on;
plot(t2, eta2);
title('Total surface elevation \eta(t), 0–300 s')
xlabel('Time (s)')
ylabel('\eta(t) (m)')
