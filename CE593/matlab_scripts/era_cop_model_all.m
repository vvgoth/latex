
clear;
g = 9.81;

%% =========================
% INPUT FILES
% =========================
filename_wind   = 'ERA5-wind-data.nc';
filename_era_hs = 'Significant height.nc';
filename_era_tm = 'era_zero_crossing_mean_period.nc';
filename_cop    = 'cop_wave_period.nc';

%% =========================
% READ WIND (ERA5)
% =========================
u10 = squeeze(ncread(filename_wind,'u10'));
v10 = squeeze(ncread(filename_wind,'v10'));

t_raw  = ncread(filename_wind,'valid_time');
time_wind = datetime(1970,1,1) + seconds(t_raw);

U = hypot(u10,v10);
windDir = mod(atan2d(-u10,-v10),360);   % FROM

%% =========================
% READ ERA5 WAVES
% =========================
Hs_era = squeeze(ncread(filename_era_hs,'shww'));
Tm_era = squeeze(ncread(filename_era_tm,'mp2'));

t_raw_era = ncread(filename_era_tm,'valid_time');
time_era  = datetime(1970,1,1) + seconds(t_raw_era);

%% =========================
% READ COPERNICUS WAVES
% =========================
Hs_cop = squeeze(ncread(filename_cop,'VHM0_WW'));
Tm_cop = squeeze(ncread(filename_cop,'VTM01_WW'));

t_raw_cop = ncread(filename_cop,'time');
time_cop  = datetime(1970,1,1) + seconds(t_raw_cop);

%% =========================
% FETCH TABLE
% =========================
Tfetch = readtable('Fetch lengths.xlsx');
fetchDir = mod(Tfetch.angle,360);
fetchLen = Tfetch.fetch_m;

fetch_hourly = zeros(size(U));
for j = 1:length(U)
    dir_bin = round(windDir(j)/7.5)*7.5;
    if dir_bin==360, dir_bin=0; end
    fetch_hourly(j) = fetchLen(fetchDir==dir_bin);
end

%% =========================
% CONTINUOUS JONSWAP HINDCAST
% =========================
N = length(U);

Hmo_model = zeros(N,1);
Tp_model  = zeros(N,1);
Hmo_star  = zeros(N,1);
t_eff     = zeros(N,1);

t_eff(1) = 3600;

[Hmo_model(1),Tp_model(1),Hmo_star(1)] = ...
    jonswap_basic(fetch_hourly(1),U(1),t_eff(1));

for j = 1:N-1

    if U(j+1)>0
        Hmo_star(j+1) = g*Hmo_model(j)/(U(j+1)^2);
    else
        Hmo_star(j+1) = 0;
    end

    t_prime_star = 368000 * Hmo_star(j+1)^(4/3);
    t_prime = t_prime_star * U(j+1) / g;

    t_eff(j+1) = t_eff(j) + t_prime;

    [Hmo_model(j+1),Tp_model(j+1)] = ...
        jonswap_basic(fetch_hourly(j+1),U(j+1),t_eff(j+1));
end

CB_coeff_Hm = 0.85; %Calibration coefficient for Hm0
CB_coeff_Tm = 0.75; %Calibration coefficient foor Tm

Hmo_model = Hmo_model * CB_coeff_Hm;
Tp_model = Tp_model * CB_coeff_Tm;
Tm_model = 0.8 * Tp_model;

%% =========================
% ALIGN ALL TIME AXES
% =========================
[time_ec, iw, ie] = intersect(time_wind, time_era);
[time_all, ic, ia] = intersect(time_ec, time_cop);

Hs_era_c = Hs_era(ie(iw(ic)));
Tm_era_c = Tm_era(ie(iw(ic)));

Hs_cop_c = Hs_cop(ia);
Tm_cop_c = Tm_cop(ia);

Hmo_c    = Hmo_model(iw(ic));
Tm_mod_c = Tm_model(iw(ic));

okH = ~isnan(Hs_era_c) & ~isnan(Hs_cop_c) & ~isnan(Hmo_c);
okT = ~isnan(Tm_era_c) & ~isnan(Tm_cop_c) & ~isnan(Tm_mod_c);

%% =========================
% STATISTICS (BIAS, RMSE, R^2)
% =========================
% --- ERA5 ---
bias_H_era = mean(Hmo_c(okH) - Hs_era_c(okH));
rmse_H_era = sqrt(mean((Hmo_c(okH) - Hs_era_c(okH)).^2));
R2_H_era   = corr(Hmo_c(okH),Hs_era_c(okH))^2;

bias_T_era = mean(Tm_mod_c(okT) - Tm_era_c(okT));
rmse_T_era = sqrt(mean((Tm_mod_c(okT) - Tm_era_c(okT)).^2));
R2_T_era   = corr(Tm_mod_c(okT),Tm_era_c(okT))^2;

% --- Copernicus ---
bias_H_cop = mean(Hmo_c(okH) - Hs_cop_c(okH));
rmse_H_cop = sqrt(mean((Hmo_c(okH) - Hs_cop_c(okH)).^2));
R2_H_cop   = corr(Hmo_c(okH),Hs_cop_c(okH))^2;

bias_T_cop = mean(Tm_mod_c(okT) - Tm_cop_c(okT));
rmse_T_cop = sqrt(mean((Tm_mod_c(okT) - Tm_cop_c(okT)).^2));
R2_T_cop   = corr(Tm_mod_c(okT),Tm_cop_c(okT))^2;

disp('--- MODEL vs ERA5 ---');
fprintf('Hs bias = %.3f m | RMSE = %.3f m | R^2 = %.3f\n', ...
        bias_H_era, rmse_H_era, R2_H_era);
fprintf('Tm bias = %.3f s | RMSE = %.3f s | R^2 = %.3f\n', ...
        bias_T_era, rmse_T_era, R2_T_era);

disp('--- MODEL vs COPERNICUS ---');
fprintf('Hs bias = %.3f m | RMSE = %.3f m | R^2 = %.3f\n', ...
        bias_H_cop, rmse_H_cop, R2_H_cop);
fprintf('Tm bias = %.3f s | RMSE = %.3f s | R^2 = %.3f\n', ...
        bias_T_cop, rmse_T_cop, R2_T_cop);

%% =========================
% PLOTS
% =========================
c_era = [0.00 0.70 0.00];
c_m   = [0.40 0.00 1.00];
c_cop = [0.83 0.03 0.03];

figure;
plot(time_all(okH),Hs_era_c(okH),'Color',c_era); hold on
plot(time_all(okH),Hs_cop_c(okH),'Color',c_cop);
plot(time_all(okH),Hmo_c(okH),'--','Color',c_m);
grid on
ylabel('H_s (m)')
title('Significant Wave Height – ERA5 vs Copernicus vs Model')
legend('ERA5','Copernicus','Model')

figure;
plot(time_all(okT),Tm_era_c(okT),'Color',c_era); hold on
plot(time_all(okT),Tm_cop_c(okT),'Color',c_cop);
plot(time_all(okT),Tm_mod_c(okT),'--','Color',c_m);
grid on
ylabel('Mean Period (s)')
title('Mean Wave Period – ERA5 vs Copernicus vs Model')
legend('ERA5','Copernicus','Model')

%% =========================
% FUNCTION
% =========================
function [H_mo, T_p, H_mo_star] = jonswap_basic(fetch,U_input,t_duration)

    g = 9.81;

    F_star = g*fetch./(U_input.^2);
    t_star = g*t_duration./U_input;

    H_mo_star_a = 0.243;
    T_p_star_a  = 8.13;

    H_mo_a = H_mo_star_a*U_input.^2/g;
    T_p_a  = T_p_star_a *U_input/g;

    F_star_eff = (t_star/68.8).^(3/2);

    H_mo_star_b = 0.0016*min(F_star,F_star_eff).^(1/2);
    T_p_star_b  = 0.286 *min(F_star,F_star_eff).^(1/3);

    H_mo_b = H_mo_star_b*U_input.^2/g;
    T_p_b  = T_p_star_b *U_input/g;

    H_mo = min(H_mo_a,H_mo_b);
    T_p  = min(T_p_a,T_p_b);
    H_mo_star = min(H_mo_star_a,H_mo_star_b);
end
