% Parameters
global g
g = 9.81;
T_all = [0.5, 0.5, 2, 2, 8, 10, 14];
h_all = [0.1, 0.5, 0.1, 0.5, 10, 50, 0.5];
Cases = (1:length(T_all))';

% Arrays
L_disp_res = zeros(length(Cases),1);
L_Guo_res = zeros(length(Cases),1);
L_Beji_res = zeros(length(Cases),1);
Guo_err = zeros(length(Cases),1);
Beji_err = zeros(length(Cases),1);

for j = 1:length(T_all)
    T = T_all(j);
    h = h_all(j); 

    % Dispersion Equation
    L_disp = dispersion_equation(T, h);
    L_disp_res(j) = L_disp;
    
    % Guo (2002)
    L_Guo = guo_method(T, h);
    L_Guo_res(j) = L_Guo;
    Guo_err(j) = error(L_Guo, L_disp);

    % Beji (2013)
    L_Beji = beji_method(T, h);
    L_Beji_res(j) = L_Beji;
    Beji_err(j) = error(L_Beji, L_disp);
end

Results = table(Cases, L_disp_res, L_Guo_res, Guo_err, L_Beji_res, Beji_err, 'VariableNames', ... 
    {'Case', 'Dispersion_L', 'Guo_L', 'Guo_Error_Percent', 'Beji_L', 'Beji_Error_Percent'});
disp(Results)

% Plotting
d_fixed = 10;
T_grid = linspace(0.1,20,400);
y = T_grid .* sqrt(g / d_fixed);
x1 = zeros(size(T_grid));
x2 = zeros(size(T_grid));
x3 = zeros(size(T_grid));

for i = 1:numel(T_grid)
    L_D = dispersion_equation(T_grid(i), d_fixed);
    L_G = guo_method(T_grid(i), d_fixed);
    L_B = beji_method(T_grid(i), d_fixed);

    x1(i) = L_D / d_fixed;
    x2(i) = L_G / d_fixed;
    x3(i) = L_B / d_fixed;
end

% Create the plot
figure; hold on; grid on; box on
plot(x1, y, 'c-', 'LineWidth', 3.5, 'DisplayName', 'Dispersion Equation');
plot(x2, y, 'y--', 'LineWidth', 2, 'DisplayName', 'Guo Method');
plot(x3, y, 'm-.', 'LineWidth', 2, 'DisplayName', 'Beji Method');

xlabel('Wavelength');
ylabel('Period');
xlim([0 16]);
ylim([0 16]);
legend('Location','northwest');


% Functions
function L = dispersion_equation(T, h)
    global g
    syms f(x)
    tolerance = 0.0001;
    err = tolerance + 1;
    L_guess = 1;
    
    while err >= tolerance 
        f_x = L_guess - ( (g / (2 * pi)) * T^2 * tanh((2 * pi / L_guess) * h ) );
        f_xt = 1 + (g / (2 * pi)) * T^2 * (2 * pi * h / L_guess^2) * sech((2 * pi / L_guess) * h)^2;
        L_new = L_guess - (f_x/f_xt);
        err = abs((L_new-L_guess)/L_new);
        L_guess = L_new;
    end
    
    L = L_guess;
end

function L = guo_method(T, h)
    global g
    b = 2.4908;
    sig = 2 * pi / T;
    
    x = (sig * h) / sqrt(g * h);
    y = x^2 * (1 - exp(-x^b))^(-1 / b);
    
    k = y / h;
    L = 2 * pi / k;
end

function L = beji_method(T, h)
    global g
    k_0 = (2 * pi / T)^2 / g;
    mu_0 = k_0 * h;
    
    fc = mu_0^1.3 * exp(-(1.1 + 2*mu_0));
    mu = mu_0 * (1 + fc) / sqrt(tanh(mu_0));
    
    k = mu / h;
    L = 2 * pi / k;
end

function error = error(L_other, L_disp)
    error = abs((L_other - L_disp) / L_disp) * 100;
end
