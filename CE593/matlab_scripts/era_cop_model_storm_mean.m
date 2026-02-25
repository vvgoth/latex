clc
clear
g = 9.81;

%% =========================
% USER SETTINGS
% =========================
mean_dur = 6;        % hours (last X-hour mean)
U_thr   = 3.0;       % storm wind threshold (m/s)
dir_thr = 45;        % direction change threshold (deg)

%% =========================
% COLORS 
% =========================
c_era = [0.00 0.70 0.00];   % ERA5
c_m   = [0.40 0.00 1.00];   % Model
c_cop = [0.83 0.03 0.03];   % Copernicus

%% =========================
% INPUT FILES
% =========================
filename_wind   = 'ERA5-wind-data.nc';
filename_era_hs = 'Significant height.nc';
filename_era_tm = 'Mean wave period.nc';
filename_cop    = 'cop_wave_period.nc';

%% =========================
% READ WIND
% =========================
u10 = squeeze(ncread(filename_wind,'u10'));
v10 = squeeze(ncread(filename_wind,'v10'));

time_wind = datetime(1970,1,1) + seconds(ncread(filename_wind,'valid_time'));

U = hypot(u10,v10);
windDir = mod(atan2d(-u10,-v10),360);

%% =========================
% READ WAVES
% =========================
Hs_era = squeeze(ncread(filename_era_hs,'shww'));
Tm_era = squeeze(ncread(filename_era_tm,'mpww'));
time_era = datetime(1970,1,1) + seconds(ncread(filename_era_tm,'valid_time'));

Hs_cop = squeeze(ncread(filename_cop,'VHM0_WW'));
Tm_cop = squeeze(ncread(filename_cop,'VTM01_WW'));
time_cop = datetime(1970,1,1) + seconds(ncread(filename_cop,'time'));

%% =========================
% FETCH TABLE
% =========================
Tfetch   = readtable('Fetch lengths.xlsx');
fetchDir = mod(Tfetch.angle,360);
fetchLen = Tfetch.fetch_m;

%% =========================
% EFFECTIVE WIND (LAST X HOURS MEAN)
% =========================
N = length(U);
U_eff   = zeros(N,1);
Dir_eff = zeros(N,1);
Fetch_eff = zeros(N,1);

for j = 1:N
    j0 = max(1,j-(mean_dur-1));

    ubar = mean(u10(j0:j));
    vbar = mean(v10(j0:j));

    U_eff(j)   = hypot(ubar,vbar);
    Dir_eff(j) = mod(atan2d(-ubar,-vbar),360);

    dir_bin = round(Dir_eff(j)/7.5)*7.5;
    if dir_bin==360, dir_bin=0; end
    Fetch_eff(j) = fetchLen(fetchDir==dir_bin);
end

%% =========================
% CONTINUOUS JONSWAP MODEL
% =========================
Hmo_model = zeros(N,1);
Tp_model  = zeros(N,1);
t_eff     = zeros(N,1);

t_eff(1) = 3600;
[Hmo_model(1),Tp_model(1)] = ...
    jonswap_basic(Fetch_eff(1),U_eff(1),t_eff(1));

for j = 1:N-1
    t_eff(j+1) = t_eff(j) + 3600;
    [Hmo_model(j+1),Tp_model(j+1)] = ...
        jonswap_basic(Fetch_eff(j+1),U_eff(j+1),t_eff(j+1));
end

CB_coeff_Hm = 0.85; %Calibration coefficient for Hm0
CB_coeff_Tm = 0.75; %Calibration coefficient foor Tm

Hmo_model = Hmo_model * CB_coeff_Hm;
Tp_model = Tp_model * CB_coeff_Tm;
Tm_model = 0.8 * Tp_model;

%% =========================
% ALIGN TIME AXES
% =========================
[time_ec, iw, ie] = intersect(time_wind,time_era);
[time_all, ic, ia] = intersect(time_ec,time_cop);

Hmo_c    = Hmo_model(iw(ic));
Tm_mod_c = Tm_model(iw(ic));

Hs_era_c = Hs_era(ie(iw(ic)));
Tm_era_c = Tm_era(ie(iw(ic)));

Hs_cop_c = Hs_cop(ia);
Tm_cop_c = Tm_cop(ia);

U_c       = U_eff(iw(ic));
windDir_c = Dir_eff(iw(ic));

%% =========================
% STORM IDENTIFICATION
% =========================
dDir = abs(diff(windDir_c));
dDir(dDir>180) = 360 - dDir(dDir>180);
dDir = [0; dDir];

isStorm = (U_c > U_thr) & (dDir < dir_thr);

%% =========================
% GROUP STORMS
% =========================
storm_id = zeros(size(isStorm));
sid = 0;
for i = 2:length(isStorm)
    if isStorm(i) && ~isStorm(i-1)
        sid = sid + 1;
    end
    if isStorm(i)
        storm_id(i) = sid;
    end
end
nStorms = sid;

%% =========================
% BUILD STORM STRUCT
% =========================
Storms = struct();
for s = 1:nStorms
    idx = storm_id == s;

    Storms(s).time = time_all(idx);

    Storms(s).Hs_model = Hmo_c(idx);
    Storms(s).Tm_model = Tm_mod_c(idx);

    Storms(s).Hs_era = Hs_era_c(idx);
    Storms(s).Tm_era = Tm_era_c(idx);

    Storms(s).Hs_cop = Hs_cop_c(idx);
    Storms(s).Tm_cop = Tm_cop_c(idx);

    Storms(s).Hs_max = max(Storms(s).Hs_model);
end

%% =========================
% GLOBAL STORM STATISTICS
% =========================
Hs_mod_all = vertcat(Storms.Hs_model);
Hs_era_all = vertcat(Storms.Hs_era);
Hs_cop_all = vertcat(Storms.Hs_cop);

Tm_mod_all = vertcat(Storms.Tm_model);
Tm_era_all = vertcat(Storms.Tm_era);
Tm_cop_all = vertcat(Storms.Tm_cop);

disp('===== GLOBAL STORM STATISTICS =====')
stats_print('Hs – ERA5',Hs_mod_all,Hs_era_all)
stats_print('Hs – Copernicus',Hs_mod_all,Hs_cop_all)
stats_print('Tm – ERA5',Tm_mod_all,Tm_era_all)
stats_print('Tm – Copernicus',Tm_mod_all,Tm_cop_all)

%% =========================
% TOP 3 STORMS – PLOTS
% =========================
[~,ord] = sort([Storms.Hs_max],'descend');
top3 = ord(1:min(3,nStorms));

for k = 1:length(top3)
    s = top3(k);

    figure
    plot(Storms(s).time,Storms(s).Hs_model,'Color',c_m,'LineWidth',2); hold on
    plot(Storms(s).time,Storms(s).Hs_era,'Color',c_era,'LineWidth',1.5)
    plot(Storms(s).time,Storms(s).Hs_cop,'Color',c_cop,'LineWidth',1.5)
    grid on
    title(sprintf('Storm %d – H_s',s))
    legend('Model','ERA5','Copernicus','Location','northwest')

    figure
    plot(Storms(s).time,Storms(s).Tm_model,'Color',c_m,'LineWidth',2); hold on
    plot(Storms(s).time,Storms(s).Tm_era,'Color',c_era,'LineWidth',1.5)
    plot(Storms(s).time,Storms(s).Tm_cop,'Color',c_cop,'LineWidth',1.5)
    grid on
    title(sprintf('Storm %d – T_m',s))
    legend('Model','ERA5','Copernicus','Location','northwest')
end

%% =========================
% FUNCTIONS
% =========================
function stats_print(name,model,data)
    ok = ~isnan(model) & ~isnan(data);

    bias = mean(model(ok)-data(ok));
    rmse = sqrt(mean((model(ok)-data(ok)).^2));

    R = corrcoef(model(ok),data(ok));
    R2 = R(1,2)^2;

    fprintf('%s | Bias %.3f | RMSE %.3f | R^2 %.3f\n', ...
        name,bias,rmse,R2);
end

function [H_mo, T_p] = jonswap_basic(fetch,U_input,t_duration)
    g = 9.81;

    F_star = g*fetch./U_input.^2;
    t_star = g*t_duration./U_input;

    H_a = 0.243*U_input.^2/g;
    T_a = 8.13*U_input/g;

    F_eff = (t_star/68.8).^(3/2);
    F_use = min(F_star,F_eff);

    H_b = 0.0016*F_use.^(1/2).*U_input.^2/g;
    T_b = 0.286 *F_use.^(1/3).*U_input/g;

    H_mo = min(H_a,H_b);
    T_p  = min(T_a,T_b);
end
