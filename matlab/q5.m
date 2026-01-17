%Load time series (1 column, 10 Hz) if available
data_path = 'ce293_2526_ts.dat';
if isfile(data_path)
    fs = 10;
    eta = readmatrix(data_path);
    eta = eta(:).';
    eta = eta(isfinite(eta));
    dt = 1 / fs;
    t = (0:numel(eta)-1) * dt;
    Nt = numel(eta);
else
    %Randomized frequencies
    N = 1000;
    edges = linspace(0.1, 2, N+1);
    widths = diff(edges);
    df = widths(1);
    f = edges(1:end-1) + widths .* rand(1, N);

    %Spectrum amplitude
    a = sqrt(2 * df .* (0.205 * 0.15^2 * 1.6^-4 .* f.^-5 .* exp(-0.75 * ...
        (1.6 .* f).^-4)));
    epsilon = 2 * pi * rand(1,N);

    %Time range
    t = 0 : 0.05 : 1200;
    Nt = length(t);

    %Calculating wave surface
    eta = sum(a.' .* cos(2 * pi * f.' .* t + epsilon.'), 1);
end

eta_mean = mean(eta);

%Mean sea level correction
if eta_mean ~= 0
    eta = eta - eta_mean;
end

%Basic statistics
eta_rms = sqrt((1/Nt) * sum(eta.^2));
skewness = ((1/eta_rms^3) * (1/Nt) * sum(eta.^3));
kurtosis = (1/eta_rms^4) * (1/Nt) * sum(eta.^4);

%Zero-up-crossing method
idx_cross = [];
tu = [];

for n_idx = 1:(Nt-1)
    if eta(n_idx) < 0 && eta(n_idx+1) > 0
        idx_cross(end+1) = n_idx;
        tu(end+1) = t(n_idx) + (-eta(n_idx) / (eta(n_idx+1) - ...
            eta(n_idx))) * (t(n_idx+1) - t(n_idx));
    end
end

n = numel(idx_cross) - 1;
Hi = zeros(1, n);
Ti = zeros(1, n);

for j = 1:n
    wave = eta(idx_cross(j):idx_cross(j+1));
    Hi(j) = max(wave) - min(wave);
    Ti(j) = tu(j+1) - tu(j);
end

[Hi_sorted, idx_sort] = sort(Hi, 'descend');
Ti_sorted = Ti(idx_sort);

%Characteristic wave parameters
Hrms = sqrt((1/n) * sum(Hi.^2));
Hm = (1/n) * sum(Hi);
Tm = (1/n) * sum(Ti);
H13 = mean(Hi_sorted(1:floor(n/3)));
T13 = mean(Ti_sorted(1:floor(n/3)));

%Histogram
x = Hi / Hm;
dx = 0.1;
edges = 0:dx:(max(x)+dx);
[counts, ~] = histcounts(x, edges);
centers = edges(1:end-1) + dx/2;

figure(1); 
bar(centers, counts, 'BarWidth', 1, ...
    'FaceColor', [0.55 0.6 1.0], ...
    'EdgeColor', [0.4 0.45 0.85]);
grid on;
xlabel('H/H_m');
ylabel('Number of Waves');
title('Non-dimensional Wave Height Histogram');

%Rayleigh comparison
Nw = numel(x);
pdf_nd = counts / (Nw * 0.1);

x_ray = linspace(0, max(x)*1.2, 200);
pdf_ray = (pi/2) .* x_ray .* exp(-pi/4 .* x_ray.^2);

figure(2);
plot(centers, pdf_nd, 'o--', 'Color', [0.5 0.7 1.0], ...
    'MarkerFaceColor', [0.5 0.7 1.0],'LineWidth', 1.2, ...
    'DisplayName','Empirical pdf');  
hold on;
plot(x_ray, pdf_ray, 'Color', [1.0 0.6 0.8],'LineWidth', 2, ...
    'DisplayName','Rayleigh pdf');
grid on;
xlabel('H/H_m');
ylabel('PDF');
title('PDF of H/H_m and Rayleigh distribution');
legend('Location','northeast');

n_rel = counts / Nw;
F_emp = cumsum(n_rel);
Q_emp = 1 - F_emp;
valid = Q_emp > 0;
x_lin = centers(valid);
sqrt_lnQ = sqrt(-log(Q_emp(valid)));

alpha = sqrt(pi)/2;
y_ray = alpha * x_lin;

p = polyfit(x_lin, sqrt_lnQ, 1);
y_fit = polyval(p, x_lin);
Rmat = corrcoef(sqrt_lnQ, y_fit);
R2 = Rmat(1,2)^2;

figure(3);
plot(x_lin, sqrt_lnQ, 'o', 'MarkerFaceColor', [0.5 0.7 1.0], ...
    'MarkerEdgeColor','none'); hold on;
plot(x_lin, y_ray, '--', 'LineWidth', 1.8, 'Color', [1.0 0.6 0.8]);                   
grid on;
xlabel('H/H_m');
ylabel('(-ln Q)^{0.5}');
title('Rayleigh Distribution Linearity Test');

text(0.1, max(sqrt_lnQ)*0.9, ...
     sprintf('y = %.4fx\n R^2 = %.4f', p(1), R2), 'Color', 'w');

%Joint distribution
tau = Ti / Tm;
tau_edges = 0:dx:(max(tau)+dx);

[joint_counts, ~, ~] = histcounts2(x, tau, edges, tau_edges);
tau_centers = tau_edges(1:end-1) + dx/2;
p_joint = joint_counts / (n * dx * dx);

Counts_TH = joint_counts.';

H_labels = arrayfun(@(a,b) sprintf('%.1f-%.1f', a, b), ...
    edges(1:end-1), edges(2:end), 'UniformOutput', false);

T_labels = arrayfun(@(a,b) sprintf('%.1f-%.1f', a, b), ...
    tau_edges(1:end-1), tau_edges(2:end), 'UniformOutput', false);

rowTotals = sum(Counts_TH, 2);
colTotals = sum(Counts_TH, 1);
grandTotal = sum(rowTotals);

JointCountsTable = array2table(Counts_TH, 'RowNames', T_labels, ...
    'VariableNames', H_labels);

JointCountsTable.("Total_H") = rowTotals;
TotalRow = array2table([colTotals, grandTotal], ...
    'RowNames', "Total_T", 'VariableNames', [H_labels, {'Total_H'}]);

JointCountsTable = [JointCountsTable; TotalRow];

figure('Name','Joint Distribution Table','NumberTitle','off');

table = uitable('Data', JointCountsTable{:,:}, ...
            'ColumnName', JointCountsTable.Properties.VariableNames, ...
            'RowName', JointCountsTable.Properties.RowNames, ...
            'Units', 'normalized', ...
            'Position',[0 0 1 1]);

table.FontSize = 12;
numCols = size(JointCountsTable{:,:},2);
table.ColumnWidth = repmat({60}, 1, numCols);

%Contour plot
figure(5);
hold on;
[T_grid, X_grid] = meshgrid(tau_centers, centers);
[~, h] = contour(T_grid, X_grid, p_joint, 5, 'LineWidth', 2);

colormap(spring);
cb = colorbar;
cb.Label.String = 'p(x,\tau)';

scatter(tau, x, 30, 'w');

xlabel('T/T_m');
ylabel('H/H_m');
title('Joint distribution of H and T');
grid on;

%Variance density spectrum
dt = t(2) - t(1);
Fs = 1/dt; %Sampling freq
Nfft = length(eta);

E = fft(eta);
E_half = E(1:floor(Nfft/2));

f_fft = (0:length(E_half)-1) * (Fs/Nfft);

S = (2/(Nfft * Fs)) * abs(E_half).^2;

m0 = trapz(f_fft, S);
Hm0 = 4 * sqrt(m0);

figure(6);
plot(f_fft, S, 'Color', [0.65 0.75 1.0], 'LineWidth', 1.8);
xlabel('Frequency (Hz)');
ylabel('S(f)  [m^2/Hz]');
title('Variance Density Spectrum');
grid on;

fprintf('H_1/3 from time domain = %.3f m\n', H13);
fprintf('H_m0 from spectrum     = %.3f m\n', Hm0);
fprintf('Difference             = %.3f m\n', abs(H13 - Hm0));

%Recreating eta from S
df_fft = f_fft(2) - f_fft(1);
a_fft = sqrt(2 * S .* df_fft);
phi_fft = 2*pi * rand(size(a_fft));

eta_new = zeros(size(t));

for i = 1:length(f_fft)
    eta_new = eta_new + a_fft(i) * cos(2*pi*f_fft(i)*t + phi_fft(i));
end

eta_new = eta_new - mean(eta_new);

%Zero-up crossing for new series
idx_cross_new = [];
tu_new = [];

for n_idx = 1:(Nt-1)
    if eta_new(n_idx) < 0 && eta_new(n_idx+1) > 0
        idx_cross_new(end+1) = n_idx;
        tu_new(end+1) = t(n_idx) + (-eta_new(n_idx) / ...
        (eta_new(n_idx+1) - eta_new(n_idx))) * (t(n_idx+1) - t(n_idx));
    end
end

n_new = numel(idx_cross_new) - 1;
Hi_new = zeros(1, n_new);
Ti_new = zeros(1, n_new);

for j = 1:n_new
    segment_new = eta_new(idx_cross_new(j):idx_cross_new(j+1));
    Hi_new(j) = max(segment_new) - min(segment_new);
    Ti_new(j) = tu_new(j+1) - tu_new(j);
end

[Hi_new_sorted, idx_sort_new] = sort(Hi_new, 'descend');
Ti_new_sorted = Ti_new(idx_sort_new);

Hrms_new = sqrt((1/n_new) * sum(Hi_new.^2));
Hm_new   = (1/n_new) * sum(Hi_new);
Tm_new   = (1/n_new) * sum(Ti_new);
H13_new  = mean(Hi_new_sorted(1:floor(n_new/3)));
T13_new  = mean(Ti_new_sorted(1:floor(n_new/3)));

fprintf('\nOriginal          |  New\n');
fprintf('H_1/3 = %.3f m   |  H_1/3 = %.3f m\n', H13, H13_new);
fprintf('H_m   = %.3f m   |  H_m   = %.3f m\n', Hm,  Hm_new);
fprintf('H_rms = %.3f m   |  H_rms = %.3f m\n', Hrms, Hrms_new);
fprintf('T_m   = %.3f s   |  T_m   = %.3f s\n', Tm,  Tm_new);
fprintf('T_1/3 = %.3f s   |  T_1/3 = %.3f s\n\n', T13, T13_new);

%Low-pass filter
f_full = (0:Nfft-1) * (Fs / Nfft);
E_filt = E;
fc = 1;

idx_filter = (f_full > fc) & (f_full < (Fs - fc));
E_filt(idx_filter) = 0;

eta_filt = real(ifft(E_filt));

Nsamp = round(20 / dt);

t20        = t(1:Nsamp);
eta20      = eta(1:Nsamp);
eta_filt20 = eta_filt(1:Nsamp);

figure(7);
plot(t20, eta20, 'Color', [0.5 0.7 1.0], 'LineWidth', 1.6); 
hold on;
plot(t20, eta_filt20, 'Color', [1.0 0.6 0.8], 'LineWidth', 1.6);
grid on;
xlabel('Time (s)');
ylabel('\eta (m)');
title('Original vs Filtered Surface Elevation');
legend('Original \eta(t)', 'Filtered \eta_{f<1Hz}(t)', 'Location','best');

%Monte-Carlo Hmax
Hs   = 8;
Ts   = 12;
Ns   = [100, 1000, 10000];
Nsim = 5000;

colors = [0.50 0.70 1.00;
          0.75 0.55 1.00;
          1.00 0.60 0.80];

figure(8); hold on;

for k = 1:length(Ns)
    
    Nw = Ns(k);
    Hmax_ratio = zeros(Nsim,1);
    
    for i = 1:Nsim
        U = rand(Nw,1);
        x = sqrt(-4/pi * log(1 - U));
        Hmax_ratio(i) = max(x);
    end
    
    bin_width = 0.05;
    edges = 0:bin_width:5;
    counts = histcounts(Hmax_ratio, edges);
    centers = edges(1:end-1) + bin_width/2;
    pdf_est = counts / (sum(counts) * bin_width);
    
    mean_val   = mean(Hmax_ratio);
    median_val = median(Hmax_ratio);
    [~, idxm]  = max(pdf_est);
    mode_val   = centers(idxm);
    
    plot(centers, pdf_est, 'LineWidth', 2, 'Color', colors(k,:), ...
        'DisplayName', sprintf(['N = %d (mean=%.2f, median=%.2f,' ...
        ' mode=%.2f)'], Nw, mean_val, median_val, mode_val));
end

xlabel('H_{max} / H_s');
ylabel('Probability density');
title('Monte-Carlo PDFs of H_{max} / H_s for Different N');
legend('Location','northwest');
grid on;
