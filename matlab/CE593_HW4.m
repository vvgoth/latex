clc; clear;

% Given Parameters
Hs0 = 4;
Ts = 8;
alpha0 = deg2rad(30);
h = 10;
g = 9.81;

%% Part (a)

L = dispersion_equation(Ts,h);
L0 = g * Ts^2 / (2*pi);
k = 2*pi/L;
C  = L/Ts;
C0 = L0/Ts;
Cg0 = C0/2;

n  = 0.5*(1 + (2*k*h)/sinh(2*k*h));
Cg = n*C;
Ks = sqrt(Cg0/Cg);

alpha10 = asin((C/C0)*sin(alpha0));
Kr = sqrt(cos(alpha0)/cos(alpha10));

Hs_10 = Hs0 * Ks * Kr;

fprintf('\n--- Part (a) Results ---\n');
fprintf('Ks         = %.3f\n', Ks);
fprintf('Kr         = %.3f\n', Kr);
fprintf('Hs at 10 m = %.3f m\n', Hs_10);

%% Part (b)

gamma = 3;
smax  = 10;

Nf = 10;
Nd = 16;

Tp = Ts / (1 - 0.132*(gamma + 0.2)^(-0.559));
fp = 1/Tp;

f = linspace(0.5*fp, 2.5*fp, 2000);

sigma = 0.07*(f <= fp) + 0.09*(f > fp);
r = exp(-(f - fp).^2 ./ (2*sigma.^2*fp^2));

alpha = 0.0624*(Hs0^2*fp^4)^0.2;
S = alpha * Hs0^2 * fp^4 .* f.^(-5) .* exp(-1.25*(fp./f).^4) .* gamma.^r;

m0 = trapz(f,S);
m0_bin = m0 / Nf;

df = f(2)-f(1);

T_rep = zeros(1,Nf);
idx0 = 1;

for i = 1:Nf
    E = 0;
    idx = idx0;
    while E < m0_bin && idx < length(f)
        E = E + S(idx)*df;
        idx = idx + 1;
    end
    f_seg = f(idx0:idx);
    S_seg = S(idx0:idx);
    f_rep = trapz(f_seg,f_seg.*S_seg)/trapz(f_seg,S_seg);
    T_rep(i) = 1/f_rep;
    idx0 = idx;
end

theta_rel = linspace(-pi/2, pi/2, Nd);
D = cos(theta_rel).^(2*smax);
D = D / sum(D); 

Kr2_sum = 0;
E_sum   = 0;

for i = 1:Nf
    T = T_rep(i);

    omega = 2*pi/T;
    C0 = g / omega;
    Ch = dispersion_equation(T,h)/T;

    for j = 1:Nd
        theta0 = alpha0 + theta_rel(j);

        % Snell's law
        sin_theta = (Ch/C0) * sin(theta0);
        if abs(sin_theta) < 1
            theta_h = asin(sin_theta);
            Kr = sqrt(cos(theta0)/cos(theta_h));
            Kr2_sum = Kr2_sum + Kr^2 * D(j);
            E_sum   = E_sum   + D(j);
        end
    end
end

Kr_irr = sqrt(Kr2_sum / E_sum);

Hs_10m = Hs0 * Ks * Kr_irr;

fprintf('\n--- Part (b) Results ---\n');
fprintf('Kr,irr     = %.3f\n',Kr_irr);
fprintf('Hs at 10 m = %.3f m\n',Hs_10m);

%% Part (c)

% Values read from Goda (2010) plots
Kr_irr_c  = 0.93;
Hs_10_c = Hs0 * Ks * Kr_irr_c;

fprintf('\n--- Part (c) Results ---\n');
fprintf('Kr         = %.3f\n', Kr_irr_c);
fprintf('Hs at 10 m = %.3f m\n', Hs_10_c);

%% Part (d)

slope = 1/30;

Kr_d = Kr_irr;

hL = h / L0;

H0p = Hs0 * Kr_d;

steepness = H0p / L0;

% Goda (1975) coefficients for H1/3 (Hs)
beta0   = 0.028 * steepness^(-0.38) * exp(20*slope^1.5);
beta1   = 0.52  * exp(4.2*slope);
betaMax = max(0.92, 0.32 * steepness^(-0.29) * exp(2.4*slope));

% Breaking-limited Hs at depth h
if hL >= 0.2
    Hs_10_d = Ks * H0p;
else
    Hs_10_d = min([ ...
        beta0*H0p + beta1*h, ...
        betaMax*H0p, ...
        Ks*H0p ...
    ]);
end

fprintf('\n--- Part (d) Results (Goda 1975 breaking) ---\n');
fprintf('Kr         = %.3f\n', Kr_d);
fprintf('Hs at 10 m = %.3f m\n', Hs_10_d);

%% Question 2
A      = 0.12;
h_vec = linspace(3,15,1000);

H_lin  = zeros(size(h_vec));
H_crit = zeros(size(h_vec));

for i = 1:length(h_vec)
    h = h_vec(i);

    % Dispersion
    L = dispersion_equation(Ts,h);
    k = 2*pi/L;
    C = L/Ts;

    % Group velocity
    n  = 0.5*(1 + (2*k*h)/sinh(2*k*h));
    Cg = n*C;

    % Shoaling
    Ks = sqrt(Cg0/Cg);

    % Refraction (Snell)
    alpha = asin((C/C0)*sin(alpha0));
    Kr = sqrt(cos(alpha0)/cos(alpha));

    % Linear transformed height
    H_lin(i) = Hs0 * Ks * Kr;

    % Incipient-breaking envelope (lecture)
    H_crit(i) = h * (A/(h/L0)) * (1 - exp(-1.5*pi*(h/L0)*(1 + 11*slope^(4/3))));
end

F = H_lin - H_crit;

% Find the FIRST sign change as you go from deep -> shallow
% (This corresponds to where shoaling curve meets the envelope)
idx = find(F(1:end-1).*F(2:end) <= 0, 1, 'first');

% Linear interpolation for the root
h1 = h_vec(idx);   
h2 = h_vec(idx+1);
F1 = F(idx);       
F2 = F(idx+1);

h_peak = h1 - F1*(h2-h1)/(F2-F1);   % root of F(h)=0

% Get H_peak by interpolating H_lin (or H_crit) at h_peak
H_peak = interp1(h_vec, H_lin, h_peak, 'linear');

fprintf('\n--- Q2 Results (Incipient Breaking) ---\n');
fprintf('h_1/3,peak = %.3f m\n', h_peak);
fprintf('H_1/3,peak = %.3f m\n', H_peak);

% Part (a)
h_goda = [1 2 3 4 6 8 10 12];

for j = 1:length(h_goda)

    h = h_goda(j);

    L = dispersion_equation(Ts, h);
    k = 2*pi/L;

    % Phase speeds
    C  = L / Ts;
    C0 = L0 / Ts;

    % --- Shoaling (group velocity) ---
    n  = 0.5*(1 + (2*k*h)/sinh(2*k*h));
    Cg = n*C;
    Cg0 = (C0)/2;

    Ks_2(j) = sqrt(Cg0/Cg);

    % --- Refraction (Snell) ---
    alpha_2(j) = asin((C/C0) * sin(alpha0));
    Kr_2(j)    = sqrt(cos(alpha0)/cos(alpha_2(j)));

    % Linear transformed height (for reference)
    H0_2(j) = Hs0 * Ks_2(j) * Kr_2(j);

    % Goda (1975) uses H0' (here taken as refraction-adjusted deepwater)
    H0_p = Hs0 * Kr_2(j);
    steepness = H0_p / L0;

    beta0   = 0.028 * steepness^(-0.38) * exp(20*slope^1.5);
    beta1   = 0.52  * exp(4.2*slope);
    betaMax = max(0.92, 0.32 * steepness^(-0.29) * exp(2.4*slope));

    hL_2 = h / L0;

    if hL_2 >= 0.2
        H_goda(j) = Ks_2(j) * H0_p;
    else
        H_goda(j) = min([ ...
            beta0*H0_p + beta1*h, ...
            betaMax*H0_p, ...
            Ks_2(j)*H0_p ...
        ]);
    end

end

fprintf('\n=== Q2(a) Results ===\n');
fprintf('  h(m)    Ks      Kr     alpha(deg)   H_lin(m)   H_Goda(m)\n');
for ii = 1:length(h_goda)
    fprintf('%5.0f   %6.3f  %6.3f    %8.2f     %8.3f   %8.3f\n', ...
        h_goda(ii), Ks_2(ii), Kr_2(ii), rad2deg(alpha_2(ii)), H0_2(ii), H_goda(ii));
end


figure; 
plot(h_goda, H0_2, 'o-', h_goda, H_goda, 's-'); grid on;
xlabel('Depth h (m)'); ylabel('H_{1/3} (m)');
legend('Linear (Ks Kr Hs0)', 'Goda (1975) breaking-limited', 'Location','best');
title('Q2: Significant wave height vs depth');

% Part (b)

% --- Deep-water parameters ---
H0 = Hs0;                    % deep-water significant wave height
L0 = g*Ts^2/(2*pi);          % deep-water wavelength
theta0 = alpha0;             % deep-water wave angle (rad)

% -------- Static set-up (mean water level rise) --------
% Goda formula including surface rollers

steep0 = H0 / L0;
lns = log(steep0);

% Normal-incidence setup (theta0 = 0)
zeta0_over_H0 = (0.0063 + 0.768*slope) ...
              - (0.0083 + 0.011*slope)*lns ...
              + (0.00372 + 0.0148*slope)*(lns^2);

zeta_theta0_0 = zeta0_over_H0 * H0;

% Oblique-incidence correction
p = 0.545 + 0.038*lns;
zeta_static = zeta_theta0_0 * (cos(theta0))^p;

% -------- Dynamic set-up (surf beat) --------
% Uses H0' (refraction-adjusted), no cosine correction

H0p = Hs0 * Kr_irr;     % use Kr for regular, Kr_irr for irregular refraction
h_shore = 0;            % shoreline

zeta_rms = (0.01 * H0p) / sqrt( (H0p/L0) * (1 + h_shore/H0p) );

fprintf('\n--- Q2(b) Wave set-up at shoreline (Goda) ---\n');
fprintf('Static set-up (mean) zeta      = %.3f m\n', zeta_static);
fprintf('Dynamic set-up (rms) zeta_rms  = %.3f m\n', zeta_rms);


x_over_H0 = h_vec ./ (slope * Hs0);

gamma_b = (A ./ (h_vec./L0)) .* ...
          (1 - exp(-1.5*pi*(h_vec./L0) .* (1 + 11*slope^(4/3))));

x_over_H0 = h_vec ./ (slope * Hs0);

gamma_lin  = H_lin  ./ h_vec;     % H/h from linear shoaling

figure; hold on; grid on;

plot(x_over_H0, gamma_lin, 'b-', 'LineWidth', 2);
plot(x_over_H0, gamma_b, 'r--', 'LineWidth', 2);

% Goda (1975) discrete points
h_goda = [1 2 3 4 6 8 10 12];
xg = h_goda ./ (slope * Hs0);
gamma_goda = H_goda ./ h_goda;

plot(xg, gamma_goda, 'ko', 'MarkerFaceColor','k');

xlabel('x / H_0');
ylabel('\gamma = H / h');

legend('Linear shoaling (no breaking)', ...
       'Incipient breaking envelope', ...
       'Goda (1975)', ...
       'Location','best');

title('Q2: Incipient breaking criterion (no breaking applied)');

