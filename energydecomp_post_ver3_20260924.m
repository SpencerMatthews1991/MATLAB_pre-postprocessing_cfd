% Compare 5 Methods of Drag Prediction (Separated and UPVD Frameworks)
clear; clc; close all;

folder_name = '.';
cv = 6;         % Control Volume to compare

% 1. Configuration, Normalization, and Styling
n_step = 1;     % Only process every n-th timestep (Set to 1 to read all data)

Vinf = 40;  % Freestream velocity (m/s)
rho = 1.1766;   % Freestream density
area = 0.14;    % reference area
q = 0.5 * rho * (Vinf)^2 * area;

f_shed = 10.22;           % True Vortex Shedding Frequency (Hz)
T_shed = 1 / f_shed;  % Shedding Period (s)
t_eval_start = 3;     % Physical time (s) where stable oscillation begins
num_cycles = 8;       % Exactly 5 full signals to display

% ------------------------------------------------------------------------
% DISSIPATION FORMULATION (no double-count):
%   phi_total = phi_molecular + phi_res
% We DROP the modeled eddy-viscosity dissipation (the old phi_vol used
% EFFECTIVE viscosity = molecular + modeled; adding phi_res to it double-counts
% the turbulent dissipation). This script now reads phi_molecular and phi_res.
%
% NEW STAR EXPORTS NEEDED (per-CV monitors over cv1..cv7, same 9-column format
% as the old phi_vol -- time, YAxisData, cv1..cv7):
%   - phi_molecular.csv   (molecular-viscosity dissipation)
%   - phi_res.csv         (resolved dissipation; set the monitor over ALL 7 cv
%                          boxes, not a single part -- that was the earlier bug)
%
% Phi_res_CV below is a FALLBACK: if phi_res.csv is still malformed, paste the
% 7 converged Volume Integral report values (cv1..cv7 order, Watts) here and the
% build-up will use them instead. Left as zeros so the script runs regardless.
% ------------------------------------------------------------------------
Phi_res_CV = [0, 0, 0, 0, 0, 0, 0];   % FALLBACK report values [W], cv1..cv7

% Exact Colors from Final_Exergy_Alex.py
C.Baseline = 'k';                     % Black
C.Sum      = [1 0.647 0];             % Orange
C.SumK     = [0.5 0 0.5];             % Purple
C.Visc     = 'r';                     % Red (aphi)
C.Mech     = [0 0.5 0];               % Dark Green (xm) matching Fluent #008000
C.Therm    = 'c';                     % Cyan (xth)
C.Vol      = 'b';                     % Blue (dXv)
C.Anablat  = 'm';                     % Magenta (anablat)

W.line = 1.6;   % Standardized line width for all plots

% 2. Helper Functions
function [time, val] = get_data(name, col, ref_time, folder, n_skip)
    % FIX 1: filenames do not carry a doubled '_vol' suffix. The caller already
    % passes names like 'phi_vol', so we append only '.csv' here.
    filename = fullfile(folder, [name, '.csv']);
    if exist(filename, 'file')
        data = readmatrix(filename);
        [~,data_column] = size(data);
        if data_column > 2
            t = data(:,1);
            v = data(:,col+1);
        else
            t = data(:,1);
            v = data(:,2);
        end

        % Apply Downsampling
        t = t(1:n_skip:end);
        v = v(1:n_skip:end);

        % FIX 2: STAR-CCM+ exports different signals at different sampling
        % rates, so each CSV can have a different number of rows. Resample
        % every signal onto the common reference grid.
        if numel(ref_time) > 2 && numel(t) ~= numel(ref_time)
            v = interp1(t, v, ref_time, 'linear', 'extrap');
            t = ref_time;
        end

        time = t;
        val  = v;
    else
        warning(['File missing: ', filename, '. Defaulting to 0.']);
        time = ref_time;
        val = zeros(size(ref_time));
    end
end

load_data = @(name, ref_time) get_data(name, cv, ref_time, folder_name, n_step);

% 3. Unified Data Loading Section
[time, drag_baseline] = load_data('drag', [0; 1]);
time_norm = (time - t_eval_start) / T_shed;
cd_baseline = drag_baseline / q;
eval_mask = (time_norm >= 0) & (time_norm <= num_cycles);
mean_baseline = mean(cd_baseline(eval_mask), 'omitnan');

% Near-Field Components
[~, drag_p]   = load_data('drag_p', time);
[~, drag_tau] = load_data('drag_tau', time);

% Volumetric Terms
% NOTE: phi_vol (EFFECTIVE viscosity = molecular + modeled) is NO LONGER used.
% Using it with phi_res double-counts the turbulent dissipation. Instead we
% load phi_molecular (molecular viscosity only); the resolved part comes from
% phi_res, and phi_total = phi_molecular + phi_res.
[~, phi_mol_vol] = load_data('phi_molecular', time);   % molecular dissipation
[~, theta_vol]  = load_data('theta_vol', time);
[~, ek_vol]     = load_data('ek_vol', time);
[~, aphi_vol]   = load_data('aphi_vol', time);
[~, anabla_vol] = load_data('anabla_vol', time);
[~, massspec]   = load_data('massspec_vol', time);

dek_dt = gradient(ek_vol, time);
dxv_dt = gradient(massspec, time);

% Mechanical / thermal boundary fluxes (plane-integrated totals)
[~, EA_total]   = load_data('ea_plane',   time);
[~, EV_total]   = load_data('ev_plane',   time);
[~, EP_total]   = load_data('ep_plane',   time);
[~, Etau_total] = load_data('Etau_plane', time);
[~, Xth_total]  = load_data('xth_plane',  time);

% NOTE: phi_res is NO LONGER loaded from CSV -- the monitor export is malformed.
% It now comes from the hard-coded Phi_res_CV vector at the top of the script.

Etp_total = EA_total + EV_total + EP_total - Etau_total;

% Momentum Terms (Manual Method) -- left at zero (not exported).
mom_drag_bot   = zeros(size(time));
mom_drag_top   = zeros(size(time));
mom_drag_left  = zeros(size(time));
mom_drag_right = zeros(size(time));

% Momentum Terms (Direct STAR-CCM+ Method) -- single plane monitor total
[~, MomentumDrag_total] = load_data('MomentumDrag_plane', time);

% 4. Component Pre-Calculations
% Method 1: Energy
% Total dissipation = MOLECULAR (time history at this CV) + RESOLVED (converged
% per-CV value). This is the non-double-counted formulation: the modeled
% eddy-viscosity part is intentionally excluded because phi_res replaces it.
phi_res_here = Phi_res_CV(cv);            % resolved, converged value at this CV [W]
phi_total = phi_mol_vol + phi_res_here;   % molecular + resolved
Cd_phi = (phi_total / Vinf) / q;
Cd_Etp = (Etp_total / Vinf) / q;
Cd_theta = (theta_vol / Vinf) / q;
Cd_dek = (dek_dt / Vinf) / q;
Cd_Energy = ((dek_dt + Etp_total + theta_vol + phi_total) / Vinf) / q;

% Method 2: Exergy
Exergy_Storage = dxv_dt + dek_dt;
Exergy_Mech_Flux = EA_total + EV_total + EP_total - Etau_total;
Exergy_Boundary_Flux = Xth_total + Exergy_Mech_Flux;
Exergy_Kinetic = Exergy_Mech_Flux + Exergy_Storage;

Cd_aphi = (aphi_vol / Vinf) / q;
Cd_anabla = (anabla_vol / Vinf) / q;
Cd_xth = (Xth_total / Vinf) / q;
Cd_xm = (Exergy_Mech_Flux / Vinf) / q;
Cd_Ex_Bound = (Exergy_Boundary_Flux / Vinf) / q;
Cd_Ex_Store = (Exergy_Storage / Vinf) / q;
Cd_Exergy = ((Exergy_Storage + Exergy_Boundary_Flux + aphi_vol + anabla_vol) / Vinf) / q;
Cd_Ex_Kinetic = (Exergy_Kinetic / Vinf) / q;

% Method 3A: Momentum Manual
Cd_mom_bot = (-1 * mom_drag_bot) / q;
Cd_mom_top = (-1 * mom_drag_top) / q;
Cd_mom_left = (-1 * mom_drag_left) / q;
Cd_mom_right = (-1 * mom_drag_right) / q;
Cd_Momentum = (-1 * (mom_drag_bot + mom_drag_top + mom_drag_left + mom_drag_right) / q );

% Method 3B: Momentum Direct
Cd_momDrag_bot   = nan(size(time));
Cd_momDrag_top   = nan(size(time));
Cd_momDrag_left  = nan(size(time));
Cd_momDrag_right = nan(size(time));
Cd_MomentumDirect = (-1 * MomentumDrag_total) / q;

% Method 4: UPVD (Energy-Based)
em_flux = EA_total + EV_total + EP_total;
total_anergy = aphi_vol + anabla_vol;
drag_residual_recov = ((em_flux - aphi_vol) ./ Vinf) - drag_p;
drag_UPVD = (total_anergy + em_flux) ./ Vinf + drag_residual_recov;

Cd_UPVD_Anergy = total_anergy ./ (Vinf * q);
Cd_UPVD_Mech = em_flux ./ (Vinf * q);
Cd_UPVD_Recov = drag_residual_recov ./ q;
Cd_UPVD_Total = drag_UPVD ./ q;

% ========================================================================
% PLOTTING SECTIONS
% ========================================================================

%% 1. Plot Method 1: Energy Balance
figure('Name', 'Method 1: Energy Balance', 'Position', [100, 100, 750, 520], 'Color', 'w');
hold on;
plot(time_norm, Cd_phi, 'Color', C.Visc, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_\Phi$');
plot(time_norm, Cd_Etp, 'Color', C.Mech, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\dot{E}_{TP}}$');
plot(time_norm, Cd_theta, 'Color', C.Therm, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_\Theta$');
plot(time_norm, Cd_dek, 'Color', C.Vol, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{dE_k/dt}$');
plot(time_norm, Cd_Energy, 'Color', C.Sum, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$\Sigma C_D$');
plot(time_norm, cd_baseline, 'Color', C.Baseline, 'LineStyle', '--', 'LineWidth', W.line, 'DisplayName', 'Near Field $C_d$');
title('Energy Balance Drag Decomposition');
xlabel('$t/T$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Energy Coefficient', 'FontSize', 12);
xlim([0, num_cycles]);
grid on; ax = gca; ax.GridAlpha = 0.3;
legend('Location', 'southoutside', 'Interpreter', 'latex', 'NumColumns', 4, 'FontSize', 10, 'Box', 'off');

%% 2. Plot Method 2: Exergy Balance
figure('Name', 'Method 2: Exergy Balance', 'Position', [200, 100, 750, 520], 'Color', 'w');
hold on;
plot(time_norm, Cd_aphi, 'Color', C.Visc, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\dot{\mathcal{A}}_\Phi}$');
plot(time_norm, Cd_xm, 'Color', C.Mech, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\dot{\mathcal{X}}_m}$');
plot(time_norm, Cd_xth, 'Color', C.Therm, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\dot{\mathcal{X}}_{th}}$');
plot(time_norm, Cd_anabla, 'Color', C.Anablat, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\mathcal{A}_{\nabla T}}$');
plot(time_norm, Cd_Ex_Store, 'Color', C.Vol, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\dot{\mathcal{X}}_{v(t)}}$');
plot(time_norm, Cd_Ex_Kinetic, 'Color', C.SumK, 'LineStyle', '--', 'LineWidth', W.line, 'DisplayName', '$\Sigma \mathcal{X}_k$');
plot(time_norm, Cd_Exergy, 'Color', C.Sum, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$\Sigma C_{\dot{\mathcal{X}}}$');
plot(time_norm, cd_baseline, 'Color', C.Baseline, 'LineStyle', '--', 'LineWidth', W.line, 'DisplayName', 'Near Field $C_d$');
title(sprintf('Far-Field Exergy Balance -- CV%d', cv));
xlabel('$t/T$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Exergy Coefficient', 'FontSize', 12);
xlim([0, num_cycles]);
grid on; ax = gca; ax.GridAlpha = 0.3; ax.XTick = 0:num_cycles;
legend('Location', 'southoutside', 'Interpreter', 'latex', 'NumColumns', 4, 'FontSize', 10, 'Box', 'off');

%% 4. Plot Method 3B: Momentum Balance (Direct STAR-CCM+)
figure('Name', 'Method 3B: Momentum Direct', 'Position', [300, 100, 750, 520], 'Color', 'w');
hold on;
plot(time_norm, Cd_MomentumDirect, 'Color', C.Sum, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$\Sigma C_D$');
plot(time_norm, cd_baseline, 'Color', C.Baseline, 'LineStyle', '--', 'LineWidth', W.line, 'DisplayName', 'Near Field $C_d$');
title('Momentum Balance (Solver Native)');
xlabel('$t/T$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Momentum Coefficient', 'FontSize', 12);
xlim([0, num_cycles]);
grid on; ax = gca; ax.GridAlpha = 0.3;
legend('Location', 'southoutside', 'Interpreter', 'latex', 'NumColumns', 4, 'FontSize', 10, 'Box', 'off');

%% 5. Plot Method 4: UPVD
figure('Name', 'Method 4: UPVD Energy-Based', 'Position', [400, 100, 750, 520], 'Color', 'w');
hold on;
plot(time_norm, Cd_UPVD_Anergy, 'Color', C.Visc, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\dot{\mathcal{A}}_{\mathrm{Total}}}$');
plot(time_norm, Cd_UPVD_Mech, 'Color', C.Mech, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\dot{E}_{m}}$');
plot(time_norm, Cd_UPVD_Recov, 'Color', C.Therm, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\Delta f, recov}$');
plot(time_norm, Cd_UPVD_Total, 'Color', C.Sum, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$\Sigma C_D$');
plot(time_norm, cd_baseline, 'Color', C.Baseline, 'LineStyle', '--', 'LineWidth', W.line, 'DisplayName', 'Near Field $C_d$');
title('UPVD Drag Decomposition');
xlabel('$t/T$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('UPVD Coefficient', 'FontSize', 12);
xlim([0, num_cycles]);
grid on; ax = gca; ax.GridAlpha = 0.3;
legend('Location', 'southoutside', 'Interpreter', 'latex', 'NumColumns', 4, 'FontSize', 10, 'Box', 'off');

%% ========================================================================
%  6. ENERGY DELTAS (BETWEEN CONTROL VOLUMES) - SEPARATED METHOD
%  ========================================================================
cv_stations = [1, 2, 3, 4, 5, 6, 7];
num_stations = length(cv_stations);

avg_XKE_CV   = zeros(1, num_stations);
avg_Aphi_CV  = zeros(1, num_stations);

for i = 1:num_stations
    current_cv = cv_stations(i);
    col_read = current_cv + 1;   % skip time + YAxisData aggregate column

    [~, temp_aphi] = get_data('aphi_vol', col_read, time, folder_name, n_step);
    [~, temp_ek]   = get_data('ek_vol', col_read, time, folder_name, n_step);
    [~, temp_mass] = get_data('massspec_vol', col_read, time, folder_name, n_step);
    temp_dek_dt = gradient(temp_ek, time);
    temp_dxv_dt = gradient(temp_mass, time);
    temp_xvt = temp_dek_dt + temp_dxv_dt;

    [~, temp_EA]   = get_data('ea_plane',   col_read, time, folder_name, n_step);
    [~, temp_EV]   = get_data('ev_plane',   col_read, time, folder_name, n_step);
    [~, temp_EP]   = get_data('ep_plane',   col_read, time, folder_name, n_step);
    [~, temp_Etau] = get_data('Etau_plane', col_read, time, folder_name, n_step);

    temp_Xm = temp_EA + temp_EV + temp_EP - temp_Etau;
    temp_XKE = temp_Xm + temp_xvt;

    temp_Cd_XKE  = (temp_XKE / Vinf) / q;
    temp_Cd_Aphi = (temp_aphi / Vinf) / q;

    avg_XKE_CV(i)  = mean(temp_Cd_XKE(eval_mask), 'omitnan');
    avg_Aphi_CV(i) = mean(temp_Cd_Aphi(eval_mask), 'omitnan');
end

avg_Aphi_CV = cumsum(avg_Aphi_CV);

dKE   = zeros(1, num_stations - 1);
dAphi = zeros(1, num_stations - 1);
label_names = cell(1, num_stations - 1);

for i = 1:(num_stations - 1)
    dKE(i) = avg_XKE_CV(i) - avg_XKE_CV(i+1);
    dAphi(i) = avg_Aphi_CV(i+1) - avg_Aphi_CV(i);
    label_names{i} = sprintf('CV%d \\rightarrow CV%d', i, i+1);
end

labels = categorical(label_names);
labels = reordercats(labels, label_names);

figure("Name", "Energy Transfer Imbalance", 'Position', [500, 100, 900, 520], 'Color', 'w');
hold on; box on;
b = bar(labels, [dKE' dAphi'], 'LineWidth', 1.0, 'BarWidth', 1);
b(1).FaceColor = [0.2 0.2 0.2];
b(2).FaceColor = [0.8 0.2 0.2];
b(1).EdgeColor = 'k';
b(2).EdgeColor = 'k';
ylabel('$\Delta C_{[]}$', 'Interpreter', 'latex', 'FontSize', 12);
title('Energy transfer between CV', 'FontWeight', 'bold');
ax = gca; ax.FontSize = 12; grid on; ax.GridAlpha = 0.3;
legend({'$\Delta \mathcal{X}_{KE}$ lost', '$\Delta \mathcal{A}_{\Phi}$ gained'}, ...
       'Interpreter', 'latex', 'FontSize', 11, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');

xtips1 = b(1).XEndPoints; xtips2 = b(2).XEndPoints;
ytips1 = b(1).YEndPoints; ytips2 = b(2).YEndPoints;
for i = 1:length(dKE)
    pct = (dAphi(i) - dKE(i)) / abs(dKE(i)) * 100;
    x_center = (xtips1(i) + xtips2(i)) / 2;
    y_max = max(ytips1(i), ytips2(i));
    if isnan(pct); label_text = 'n/a'; else; label_text = sprintf('%+.1f%%', pct); end
    text(x_center, y_max, label_text, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom', 'FontSize', 11, 'FontWeight', 'bold');
end
y_top = max([ytips1(:); ytips2(:)], [], 'omitnan');
y_bot = min([ytips1(:); ytips2(:); 0], [], 'omitnan');
if ~isfinite(y_top) || y_top <= y_bot; y_top = y_bot + 1; end
ylim([y_bot, y_top * 1.15]);

%% ========================================================================
%  7. FAR-FIELD BUILD-UP ALONG PLANES (time-averaged)
%  ========================================================================
avg_t_start = 3.0;
avg_t_end   = 3.45;

EA_p   = planes_meanwin('ea_plane.csv',   avg_t_start, avg_t_end);
EP_p   = planes_meanwin('ep_plane.csv',   avg_t_start, avg_t_end);
EV_p   = planes_meanwin('ev_plane.csv',   avg_t_start, avg_t_end);
EW_p = planes_meanwin('ew_plane.csv', avg_t_start, avg_t_end);
Etau_p = planes_meanwin('Etau_plane.csv', avg_t_start, avg_t_end);
EtauRes_p = planes_meanwin('Etau_res_flux.csv', avg_t_start, avg_t_end);       %flux term
Phi_mol_p = planes_meanwin('phi_molecular.csv', avg_t_start, avg_t_end);  % molecular
Phi_res_p = planes_meanwin('phi_res.csv',       avg_t_start, avg_t_end);  % resolved

% Align mechanical/molecular terms to the common number of stations.
nP = min([numel(EA_p), numel(EP_p), numel(EV_p), numel(EW_p), ...
          numel(Etau_p), numel(EtauRes_p), numel(Phi_mol_p)]);
EA_p=EA_p(1:nP); EP_p=EP_p(1:nP); EV_p=EV_p(1:nP); EW_p=EW_p(1:nP);
Etau_p=Etau_p(1:nP); EtauRes_p=EtauRes_p(1:nP);
Phi_mol_p = Phi_mol_p(1:nP);

% Resolved dissipation: prefer the per-CV phi_res.csv (proper 9-column monitor).
% If that export is still malformed (returns the wrong length), fall back to the
% hard-coded Phi_res_CV vector of converged report values.
if numel(Phi_res_p) == nP
    Phi_res_p = Phi_res_p(1:nP);
elseif numel(Phi_res_CV) >= nP
    warning('phi_res.csv gave %d value(s), not %d per-CV; using hard-coded Phi_res_CV.', numel(Phi_res_p), nP);
    Phi_res_p = Phi_res_CV(1:nP);
else
    warning('No valid phi_res source; setting resolved dissipation to 0.');
    Phi_res_p = zeros(1, nP);
    EtauRes_p = EtauRes_p(1:nP);
end

% Total dissipation per slab = MOLECULAR + RESOLVED (no modeled term -> no
% double count), THEN accumulate so Phi_n is the running dissipation to plane n.
Phi_p = Phi_mol_p + Phi_res_p;
Phi_p = cumsum(Phi_p);

% Normalize: (W / Vinf) / q -> Cd units.
nrm    = 1 / (Vinf * q);
Ek_n   = (EA_p + EV_p + EW_p) * nrm;     % longitudinal + both transverse (Eq 19)
Ep_n   = EP_p * nrm;
Etau_n = (-Etau_p + EtauRes_p) * nrm;    % mean (Eq 21) + resolved (Eq A9) shear flux
Phi_n  = Phi_p * nrm;
Total_n = Phi_n + Ek_n + Ep_n + Etau_n;

% Reference C_D over the SAME seconds window as the build-up curves
ff_mask = (time >= avg_t_start) & (time <= avg_t_end);
CD_ff   = mean(cd_baseline(ff_mask), 'omitnan');
pct_ff = 100 * mean(Total_n, 'omitnan') / CD_ff;

xf  = linspace(1, nP, 400);
pcf = @(y) pchip(pidx, y, xf);

figure('Name', 'Far-Field Build-up Along Planes', 'Position', [600, 100, 900, 540], 'Color', 'w');
hold on; box on; grid on; ax = gca; ax.GridAlpha = 0.3;
plot(xf, pcf(Ek_n),    '-', 'Color', [0 0 0.55],    'LineWidth', W.line,     'DisplayName', '$\dot{E}_{k,n}$');
plot(xf, pcf(Ep_n),    '-', 'Color', [0 0.6 0],     'LineWidth', W.line,     'DisplayName', '$\dot{E}_{p,n}$');
plot(xf, pcf(Etau_n),  '-', 'Color', [0.1 0.75 0.9],'LineWidth', W.line,     'DisplayName', '$\dot{E}_{\tau,n}$');
plot(xf, pcf(Phi_n),   '-', 'Color', [0.5 0.5 0.5], 'LineWidth', W.line,     'DisplayName', '$\phi_n$');
plot(xf, pcf(Total_n), '-', 'Color', [0.85 0 0],    'LineWidth', W.line+0.4, 'DisplayName', '$(\Phi_n + \Sigma\dot{E}_n)$');
plot([1 nP], [CD_ff CD_ff], 'k--', 'LineWidth', W.line, 'DisplayName', '$C_D$');
xlim([1 nP]);
xlabel('Plane', 'FontSize', 12);
ylabel('Normalized quantities', 'FontSize', 12);
title('IDDES: Far-field normalized quantities trend along planes (PCHIP interpolation)');
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 11, 'Box', 'on');
annotation('textbox', [0.30 0.82 0.42 0.06], ...
    'String', sprintf('Total C_D reproduced is on average %.2f%% of C_D', pct_ff), ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'FitBoxToText', 'on', 'FontSize', 10);

fprintf('Far-field build-up: mean reconstructed C_D = %.4f  (%.2f%% of near-field %.4f)\n', ...
        mean(Total_n, 'omitnan'), pct_ff, CD_ff);
fprintf('Resolved dissipation per CV (Phi_res_CV): [%s] W\n', num2str(Phi_res_CV, '%.3g '));

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================
function export_all_plots()
    figHandles = findall(0, 'Type', 'figure');
    if isempty(figHandles); disp('No figures to export.'); return; end
    export_dir = 'exported_plots';
    if ~exist(export_dir, 'dir'); mkdir(export_dir); end
    disp('Exporting plots at 600 DPI. Please wait...');
    for i = 1:length(figHandles)
        fig = figHandles(i);
        fig_name = fig.Name;
        if isempty(fig_name); fig_name = sprintf('Figure_%d', fig.Number); end
        safe_name = regexprep(fig_name, '[\\\/\:\*\?\"\<\>\|]', '_');
        safe_name = strrep(safe_name, ' ', '_');
        filepath = fullfile(export_dir, [safe_name, '.png']);
        exportgraphics(fig, filepath, 'Resolution', 600);
        fprintf('Saved: %s\n', filepath);
    end
    disp('All plots exported successfully!');
end

function v = planes_meanwin(fname, t_start, t_end)
    % Time-average each plane/CV column over [t_start, t_end] in SECONDS.
    % Layout: col 1 = time, col 2 = YAxisData (skipped), cols 3..end = plane1..N.
    if ~exist(fname, 'file')
        warning(['File missing: ', fname, '. Term set to 0.']);
        v = zeros(1, 7); return;
    end
    d = readmatrix(fname);
    t = d(:,1);
    if size(d,2) >= 3
        M = d(:,3:end);
    else
        M = d(:,2);
    end
    m = (t >= t_start) & (t <= t_end);
    if ~any(m)
        warning('Window [%g, %g] s has no samples in %s (data %g..%g s). Averaging all.', ...
                 t_start, t_end, fname, min(t), max(t));
        m = true(size(t));
    end
    v = mean(M(m,:), 1, 'omitnan');
end