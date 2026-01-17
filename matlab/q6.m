clc; clear;

% Parameters
N = 8;

Hs = [7.19; 9.32; 5.72; 8.11; 6.15; 6.03; 7.06; 6.37];
Ts = [11.32; 12.71; 9.95; 11.7; 10.32; 10.22; 10.77; 10.51];

Hs_desc = sort(Hs, 'descend');

i = (1:N)';

% Gumbel FT-Ia
alpha_g = 0;
beta_g  = 1;

F_g = 1 - (i - alpha_g) ./ (N + beta_g);
y_g = -log( -log(F_g) );

x = Hs_desc;

Ybar_g = (1/N) * sum(y_g);
Xbar_g = (1/N) * sum(x);

VarY_g = (1/N) * sum( (y_g - Ybar_g).^2 );
VarX_g  = (1/N)*sum( (x - Xbar_g).^2 );
CovYX_g = (1/N) * sum( (y_g - Ybar_g) .* (x - Xbar_g) );

A_g = CovYX_g / VarY_g;
B_g = Xbar_g - A_g * Ybar_g;

% Weibull
k = 2;

alpha_w = 0.20 + 0.27 / sqrt(k);
beta_w  = 0.20 + 0.23 / sqrt(k);

F_w = 1 - (i - alpha_w) ./ (N + beta_w);
y_w = ( -log(1 - F_w) ).^(1/k);

Ybar_w = (1/N) * sum(y_w);
Xbar_w = (1/N) * sum(x);

VarY_w  = (1/N) * sum( (y_w - Ybar_w).^2 );
VarX_w  = (1/N)*sum( (x - Xbar_w).^2 );
CovYX_w = (1/N) * sum( (y_w - Ybar_w) .* (x - Xbar_w) );

A_w = CovYX_w / VarY_w;
B_w = Xbar_w - A_w * Ybar_w;

%% Kolmogorov–Smirnov Test
Hs_asc = sort(Hs);
i_asc = (1:N)';

S1 = (i_asc - 1) / N;
S2 = i_asc / N;

% Gumbel FT-Ia
F_g_theo = exp( -exp( -(Hs_asc - B_g) / A_g ) );

Dg_1 = F_g_theo - S1;
Dg_2 = S2 - F_g_theo;

D_KS_g = max( max(Dg_1), max(Dg_2) );

% Weibull
F_w_theo = 1 - exp( -((Hs_asc - B_w) ./ A_w).^k );

Dw_1 = F_w_theo - S1;
Dw_2 = S2 - F_w_theo;

D_KS_w = max( max(Dw_1), max(Dw_2) );

d_crit = 0.24;

pass_g = D_KS_g < d_crit;
pass_w = D_KS_w < d_crit;

figure; hold on;

plot(Hs_asc, F_w_theo, 'b-', 'LineWidth',1.5);

for j = 1:N
    plot([Hs_asc(j) Hs_asc(j)], [S1(j) S2(j)], 'k--', 'LineWidth',1.2);
    
    if j < N
        plot([Hs_asc(j) Hs_asc(j+1)], [S2(j) S2(j)], 'k--', 'LineWidth',1.2);
    end
end

grid on;
xlabel('H_s (m)');
ylabel('Cumulative probability');
title('KS Test – Weibull (k = 1.4)');

legend('F(x) theoretical','S(x) empirical','Location','southeast');

%% R² Test
% Gumbel FT-Ia
r_g  = CovYX_g / sqrt(VarX_g*VarY_g);
R2_g = r_g^2;

% Weibull
r_w  = CovYX_w / sqrt(VarX_w*VarY_w);
R2_w = r_w^2;

%% MIR Test
v = 1;

% Gumbel FT-Ia 
a_gM = -2.364 + 0.54*(v^(5/2));
b_gM = -0.2665 - 0.0457*(v^(5/2));
c_gM = -0.044;

Delta_r_g = 1 - r_g;
Delta_rmean_g = exp(a_gM + b_gM*log(N) + c_gM*(log(N))^2);

MIR_g = Delta_r_g / Delta_rmean_g;

% Weibull
a_wM = -2.160 + 0.113*v;
b_wM = -0.3788 - 0.0979*v;
c_wM = -0.041;

Delta_r_w = 1 - r_w;
Delta_rmean_w = exp(a_wM + b_wM*log(N) + c_wM*(log(N))^2);

MIR_w = Delta_r_w / Delta_rmean_w;

%% DOL Test
x1 = Hs_desc(1);
xbar = (1/N)*sum(Hs_desc);
s2 = (1/N)*sum((Hs_desc - xbar).^2);

xi = (x1 - xbar) / s2;

% Gumbel FT-Ia
% Upper bound
aU_g = -0.579 + 0.468*v;
bU_g =  1.496 - 0.227*v^2;
cU_g = -0.038;

xi95_g = aU_g + bU_g*log(N) + cU_g*(log(N))^2;

% Lower bound
aL_g = 0.257 + 0.133*v^2;
bL_g = 0.452 - 0.118*v^2;
cL_g = 0.032;

xi5_g = aL_g + bL_g*log(N) + cL_g*(log(N))^2;

pass_DOL_g = (xi > xi5_g) && (xi < xi95_g);

% Weibull
% Upper bound
aU_w = -0.322 + 0.641*sqrt(v);
bU_w =  1.414 - 0.326*v;
cU_w = -0.069;

xi95_w = aU_w + bU_w*log(N) + cU_w*(log(N))^2;

% Lower bound
aL_w = 0.050 + 0.182*sqrt(v^(3/2));
bL_w = 0.592 - 0.139*sqrt(v^(3/2));
cL_w = 0;

xi5_w = aL_w + bL_w*log(N) + cL_w*(log(N))^2;

pass_DOL_w = (xi > xi5_w) && (xi < xi95_w);

%% REC Test
% Gumbel FT-Ia
a_g = -1.444;
b_g = -0.2733 + 0.0414*sqrt(v^5);
c_g = -0.045;

Delta_r95_g = exp( a_g + b_g*log(N) + c_g*(log(N))^2 );

pass_REC_g = Delta_r_g < Delta_r95_g;

% Weibull
a_w = -1.188 + 0.073*sqrt(v);
b_w = -0.4401 - 0.0846*(v^(3/2));
c_w = -0.039;

Delta_r95_w = exp( a_w + b_w*log(N) + c_w*(log(N))^2 );

pass_REC_w = Delta_r_w < Delta_r95_w;

%% Test Results

fprintf('\n---- GOODNESS-OF-FIT RESULTS ----\n\n');

fprintf('--- Kolmogorov–Smirnov Test ---\n');
fprintf('Gumbel FT-Ia : D = %.4f  | Pass = %d\n', D_KS_g, pass_g);
fprintf('Weibull k=1.4: D = %.4f  | Pass = %d\n\n', D_KS_w, pass_w);

fprintf('--- Coefficient of Determination ---\n');
fprintf('Gumbel FT-Ia : R^2 = %.4f\n', R2_g);
fprintf('Weibull k=1.4: R^2 = %.4f\n\n', R2_w);

fprintf('--- MIR Criterion ---\n');
fprintf('Gumbel FT-Ia : MIR = %.4f\n', MIR_g);
fprintf('Weibull k=1.4: MIR = %.4f\n\n', MIR_w);

fprintf('--- DOL Criterion ---\n');
fprintf('Gumbel FT-Ia : Pass = %d\n', pass_DOL_g);
fprintf('Weibull k=1.4: Pass = %d\n\n', pass_DOL_w);

fprintf('--- REC Criterion ---\n');
fprintf('Gumbel FT-Ia : Pass = %d\n', pass_REC_g);
fprintf('Weibull k=1.4: Pass = %d\n\n', pass_REC_w);

%% Confidence bounds on probability plots
sigma_x = sqrt( (1/N)*sum((Hs_desc - mean(Hs_desc)).^2) );

z70 = norminv((1 + 0.70) / 2);

% Gumbel FT-I
y_plot_g = linspace(min(y_g), max(y_g), 200);
x_fit_g  = A_g*y_plot_g + B_g;

a1_g = 0.64;
a2_g = 9.0;
kappa_g = 0.93;
c_g = 0;
alpha_g = 1.33;

a_g = a1_g * exp( a2_g*N^(-1.3) + kappa_g*(-log(v))^2 );

sigma_z_g = sqrt( 1 + a_g*(y_plot_g - c_g + alpha_g*log(v)).^2 ) / sqrt(N);
sigma_X_g = sigma_z_g * sigma_x;

x70U_g = x_fit_g + z70*sigma_X_g;
x70L_g = x_fit_g - z70*sigma_X_g;


figure; hold on; box on
plot(y_g, Hs_desc, 'ko', 'MarkerFaceColor','k')
plot(y_plot_g, x_fit_g, 'r-', 'LineWidth',1.8)
plot(y_plot_g, x70U_g, 'k--', y_plot_g, x70L_g, 'k--')

grid on
xlabel('y = -ln[-ln(P)]')
ylabel('H_s (m)')
title('Gumbel FT-Ia Probability Plot with 70%Confidence Bounds')

legend('Data','Fitted line','70% CI','', ...
       'Location','northwest')
% Weibull
y_plot_w = linspace(min(y_w), max(y_w), 200);
x_fit_w  = A_w*y_plot_w + B_w;

a1_w = 2.05;
a2_w = 11.4;
kappa_w = 0.69;
c_w = 0.4;
alpha_w = 0.72;

a_w = a1_w * exp( a2_w*N^(-1.3) + kappa_w*(-log(v))^2 );

sigma_z_w = sqrt( 1 + a_w*(y_plot_w - c_w + alpha_w*log(v)).^2 ) / sqrt(N);
sigma_X_w = sigma_z_w * sigma_x;

x70U_w = x_fit_w + z70*sigma_X_w;
x70L_w = x_fit_w - z70*sigma_X_w;


figure; hold on; box on
plot(y_w, Hs_desc, 'ko', 'MarkerFaceColor','k')
plot(y_plot_w, x_fit_w, 'b-', 'LineWidth',1.8)
plot(y_plot_w, x70U_w, 'k--', y_plot_w, x70L_w, 'k--')

grid on
xlabel('y = [-ln(1 - F)]^{1/k}')
ylabel('H_s (m)')
title('Weibull Probability Plot with 70% Confidence Bounds')

legend('Data','Fitted line','70% CI','',  ...
       'Location','northwest')

%% 100-year return value with confidence intervals
R = 100;
F_R = 1 - 1/R;

%Gumbel FT-Ia
yR_g = -log( -log(F_R) );
XR_g = A_g*yR_g + B_g;

sigma_zR_g = sqrt( 1 + a_g*(yR_g - c_g + alpha_g*log(v))^2 ) / sqrt(N);
sigma_XR_g = sigma_zR_g * sigma_x;

XR_g_70 = [XR_g - z70*sigma_XR_g , XR_g + z70*sigma_XR_g];
XR_g_90 = [XR_g - z90*sigma_XR_g , XR_g + z90*sigma_XR_g];

% Weibull
yR_w = ( -log(1 - F_R) ).^(1/k);
XR_w = A_w*yR_w + B_w;

sigma_zR_w = sqrt( 1 + a_w*(yR_w - c_w + alpha_w*log(v))^2 ) / sqrt(N);
sigma_XR_w = sigma_zR_w * sigma_x;

XR_w_70 = [XR_w - z70*sigma_XR_w , XR_w + z70*sigma_XR_w];
XR_w_90 = [XR_w - z90*sigma_XR_w , XR_w + z90*sigma_XR_w];

% Results

fprintf('\n---- 100-YEAR RETURN VALUE ----\n\n');

fprintf('Gumbel FT-Ia:\n');
fprintf('H_s,100 = %.3f m\n', XR_g);
fprintf('70%% CI  = [%.3f , %.3f] m\n', XR_g_70(1), XR_g_70(2));
fprintf('90%% CI  = [%.3f , %.3f] m\n\n', XR_g_90(1), XR_g_90(2));

fprintf('Weibull (k = 1.4):\n');
fprintf('H_s,100 = %.3f m\n', XR_w);
fprintf('70%% CI  = [%.3f , %.3f] m\n', XR_w_70(1), XR_w_70(2));
fprintf('90%% CI  = [%.3f , %.3f] m\n\n', XR_w_90(1), XR_w_90(2));
