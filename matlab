%% rover_mass_sweep.m
% Heatmaps of allowable rover mass vs. wheel Diameter (D) and Width (b)
% Model: Bekker (pressure–sinkage) + Janosi–Hanamoto (shear/traction)
% Terrain: Moon, slope = 15 deg, 6 driven wheels
% Result: m_max = min( m_sink(z=alpha*r), m_trac(z=alpha*r, i) )

clear; clc;

%% ===================== USER SETTINGS =====================
% Wheel size ranges (edit as needed)
D_range = linspace(0.50, 1.20, 71);   % Wheel diameter [m]  (min≈Lunokhod-class, max≈Lunar Cruiser-class)
b_range = linspace(0.15, 0.50, 71);   % Wheel width [m]

% Design sinks at fixed sinkage fraction(s): z = alpha * r
alpha_list = [0.05, 0.07, 0.10];      % try 5%, 7%, 10% of radius

% Vehicle & terrain
Nw     = 6;                 % number of driven wheels
theta  = deg2rad(15);       % slope [rad]
g_moon = 1.62;              % lunar gravity [m/s^2]

% ---------------- Soil parameters (White Paper) ----------------
% Bekker pressure–sinkage: p = (kc/b + kphi) * z^n
n     = 1.0;                % Bekker exponent
kc    = 1.4e3;              % cohesive modulus [N/m^2]
kphi  = 8.3e5;              % frictional modulus [N/m^3]

% Janosi–Hanamoto / Mohr–Coulomb shear:
c     = 0.17e3;             % cohesion [Pa] = 0.17 kPa
phi   = deg2rad(38);        % internal friction angle [rad]
Ksh   = 0.015;              % shear deformation modulus [m]  (Janosi parameter)

% Mobility / margins
islip  = 0.20;              % design slip (20%)
gammaF = 1.30;              % traction safety factor on required force

% Plot style
cmapName = 'turbo';         
%% =========================================================

% Mesh
[D_mesh, b_mesh] = meshgrid(D_range, b_range);
r_mesh = 0.5 * D_mesh;

% Pre-allocate outputs for each alpha
results = struct();

for aidx = 1:numel(alpha_list)
    alpha = alpha_list(aidx);
    z_mesh = alpha .* r_mesh;   % fixed sinkage at design point
    % === Exact chord contact length (rigid wheel) ===
    Lc = sqrt(2 .* r_mesh .* z_mesh); 

    % Pressure at fixed sinkage
    p_mesh = (kc ./ b_mesh + kphi) .* (z_mesh .^ n);     % [Pa]
    
    %% --- Sinkage-limited mass ---
    % Per-wheel allowable normal load Ww_max = (kc + b*kphi) * z^n * Lc
    Ww_max = (kc + b_mesh .* kphi) .* (z_mesh .^ n) .* Lc;   % [N]
    m_sink = (Nw .* Ww_max) ./ g_moon;                       % [kg]
    
    %% --- Traction-limited mass ---
    % Shear stress tau = (c + sigma*tan(phi)) * (1 - exp(-j/K))
    sigma  = p_mesh;                     % normal stress ~ pressure
    j_mesh = islip .* Lc;                % shear displacement
    tau    = (c + sigma .* tan(phi)) .* (1 - exp(-j_mesh ./ Ksh));  % [Pa]
    tau    = max(tau, 0);                % safeguard
    
    % Available thrust (all wheels): T_avail = Nw * tau * b * Lc
    T_avail = Nw .* tau .* b_mesh .* Lc; % [N]
    
    % === Rolling resistance rising with sinkage (empirical) ===
    z_over_r = z_mesh ./ r_mesh;           % sinkage fraction
    % Crr(z/r) = base + a1*(z/r) + a2*(z/r)^2 (calibrated for soft regolith)
    Crr_mesh = 0.05 + 1.0 .* z_over_r + 2.0 .* (z_over_r.^2);
    Crr_mesh = min(Crr_mesh, 0.25);        % saturate at 0.25

    % Required force up slope with rolling: F_req = m g (sinθ + Crr cosθ)
    denom = gammaF * g_moon .* (sin(theta) + Crr_mesh .* cos(theta)); % [N/kg]
    m_trac = T_avail ./ denom;                                  % [kg]
    
    % Governing allowable mass
    m_max = min(m_sink, m_trac);
    
    % Save
    results(aidx).alpha = alpha;
    results(aidx).m_sink = m_sink;
    results(aidx).m_trac = m_trac;
    results(aidx).m_max  = m_max;
    results(aidx).D_mesh = D_mesh;
    results(aidx).b_mesh = b_mesh;
    results(aidx).Crr    = Crr_mesh;    
end

%% ======== Plot: m_max heatmap for each alpha, with constraint boundary ========
for aidx = 1:numel(alpha_list)
    alpha = results(aidx).alpha;
    m_sink = results(aidx).m_sink;
    m_trac = results(aidx).m_trac;
    m_max  = results(aidx).m_max;
    Dm     = results(aidx).D_mesh;
    bm     = results(aidx).b_mesh;

    figure('Color','w');
    % Heatmap
    contourf(Dm, bm, m_max, 24, 'LineColor', 'none');
    colormap(cmapName); colorbar; hold on;
    xlabel('Wheel Diameter D [m]'); ylabel('Wheel Width b [m]');
    title(sprintf('Allowable Rover Mass m_{max} [kg] at \\alpha = %.2f (z = \\alpha r)', alpha));
    set(gca,'FontSize',11); grid on; box on;

    % Constraint boundary where m_sink = m_trac
    % hold on;
    % C = contour(Dm, bm, m_sink - m_trac, [0 0], 'k', 'LineWidth', 1.5);
    % if ~isempty(C)
    %     % Add a text label manually near the midpoint of the contour
    %     clabel(C, 'Color','k', 'FontSize',9, 'LabelSpacing', 300);
    %     % Optional manual annotation
    %     midIdx = round(size(C,2)/2);
    %     text(C(1,midIdx), C(2,midIdx), 'm_{sink}=m_{trac}', ...
    %         'Color','k', 'FontSize',9, 'FontWeight','bold', ...
    %         'BackgroundColor','w', 'Margin',1);
    % end

    % Mark the optimum point
    [mmax_val, idx] = max(m_max(:));
    [i_opt, j_opt]  = ind2sub(size(m_max), idx);
    D_opt = Dm(i_opt, j_opt); b_opt = bm(i_opt, j_opt);
    plot(D_opt, b_opt, 'wo', 'MarkerFaceColor','k', 'MarkerSize', 7);
    text(D_opt, b_opt, sprintf('  Max = %.0f kg', mmax_val), 'Color','k', 'FontSize',10, 'FontWeight','bold');

    % Annotate governing region
    text(min(D_range)+0.02, max(b_range)-0.02, 'Left of curve: traction governs', 'Color',[0 0 0], 'FontSize',9, 'VerticalAlignment','top');
    text(min(D_range)+0.02, max(b_range)-0.06, 'Right of curve: sinkage governs',  'Color',[0 0 0], 'FontSize',9, 'VerticalAlignment','top');
end

%% ======== Optional: show separate m_sink and m_trac maps for one alpha ========
show_detail_for_alpha = alpha_list(min(2,numel(alpha_list))); % pick the 2nd alpha by default, if present
ai = find(abs([results.alpha] - show_detail_for_alpha) < 1e-9, 1);

if ~isempty(ai)
    Dm = results(ai).D_mesh; bm = results(ai).b_mesh;
    figure('Color','w');
    subplot(1,2,1);
    contourf(Dm, bm, results(ai).m_sink, 24, 'LineColor','none'); colormap(cmapName); colorbar;
    xlabel('D [m]'); ylabel('b [m]'); title(sprintf('m_{sink} [kg] (\\alpha=%.2f)', results(ai).alpha));
    set(gca,'FontSize',11); grid on; box on;

    subplot(1,2,2);
    contourf(Dm, bm, results(ai).m_trac, 24, 'LineColor','none'); colormap(cmapName); colorbar;
    xlabel('D [m]'); ylabel('b [m]'); title(sprintf('m_{trac} [kg] (\\alpha=%.2f)', results(ai).alpha));
    set(gca,'FontSize',11); grid on; box on;
end

%% ======== Console summary of the best design(s) ========
fprintf('==== Summary (per alpha) ====\n');
for aidx = 1:numel(alpha_list)
    alpha = results(aidx).alpha;
    m_max  = results(aidx).m_max;
    Dm     = results(aidx).D_mesh;
    bm     = results(aidx).b_mesh;
    [mmax_val, idx] = max(m_max(:));
    [i_opt, j_opt]  = ind2sub(size(m_max), idx);
    fprintf('alpha = %.2f: Max allowable mass = %.0f kg at D = %.2f m, b = %.2f m\n', ...
        alpha, mmax_val, Dm(i_opt,j_opt), bm(i_opt,j_opt));
end

%% ======== Notes ========
% - m_sink uses Ww_max = (kc + b*kphi) * z^n * Lc, with z = alpha*r and Lc = sqrt(2*r*z).
% - m_trac uses tau = (c + p*tan(phi)) * (1 - exp(-(i*Lc)/Ksh)), p = (kc/b + kphi) * z^n,
%   T_avail = Nw * tau * b * Lc, and F_req = gammaF * m g (sin(theta) + Crr cos(theta)).
% - m_max = min(m_sink, m_trac).
% - Increase b or D generally increases m_max; more wheels (Nw) also helps.
% - Re-run with softer soil (lower c, phi, Ksh; higher Crr) to stress-test design.

%% ======== Sinkage-dependent rolling resistance + torque demand ========
% Model: Crr(z/r) = Crr_min + k1*(z/r) + k2*(z/r).^2   (empirical, tunable)
% Rationale: captures the sharp rise in rolling/compaction/bulldozing losses with sinkage.
% Calibrate k1,k2 to your data; defaults below are conservative for fluffy regolith.

Crr_min = 0.05;   % base rolling on firmed soil (dimensionless)
k1      = 0.60;   % linear rise with z/r (dimensionless)
k2      = 0.00;   % optional curvature term (start at 0)

for aidx = 1:numel(alpha_list)
    alpha = results(aidx).alpha;

    Dm = results(aidx).D_mesh;
    bm = results(aidx).b_mesh;
    r_mesh = 0.5 .* Dm;

    % sinkage at this alpha and contact length already computed earlier:
    z_mesh = alpha .* r_mesh;
    Lc     = sqrt(2 .* r_mesh .* z_mesh); %#ok<NASGU>  % (kept for reference)

    % --- rolling resistance coefficient as function of sinkage fraction ---
    z_over_r = z_mesh ./ r_mesh;     % = alpha, but keep explicit for clarity
    Crr_z = Crr_min + k1 .* z_over_r + k2 .* (z_over_r.^2);
    Crr_z = max(Crr_z, Crr_min);     % never below base

    % Store map
    results(aidx).Crr_z = Crr_z;

    % --- torque demand per wheel at governing mass m_max ---
    m_max = results(aidx).m_max;

    % Along-slope components
    F_grade = m_max .* g_moon .* sin(theta);                 % [N]
    F_roll  = m_max .* g_moon .* (Crr_z .* cos(theta));      % [N]
    F_total = F_grade + F_roll;                              % [N]

    % Per-wheel torque at tire (no gear/efficiency yet)
    T_wheel = (r_mesh .* F_total) ./ Nw;                     % [N·m]

    results(aidx).T_wheel = T_wheel;

    % ---- PLOTS for this alpha ----
    figure('Color','w');
    subplot(1,2,1);
    contourf(Dm, bm, Crr_z, 24, 'LineColor','none'); colormap(turbo); cb=colorbar;
    cb.Label.String = 'C_{rr}(z/r)'; grid on; box on;
    xlabel('Wheel Diameter D [m]'); ylabel('Wheel Width b [m]');
    title(sprintf('Rolling Resistance Coefficient C_{rr} at \\alpha = %.2f', alpha));
    set(gca,'FontSize',11);

    subplot(1,2,2);
    contourf(Dm, bm, T_wheel, 24, 'LineColor','none'); colormap(turbo); cb=colorbar;
    cb.Label.String = 'Per-wheel torque T_{wheel} [N·m]'; grid on; box on;
    xlabel('Wheel Diameter D [m]'); ylabel('Wheel Width b [m]');
    title(sprintf('Wheel Torque Demand at m_{max} (\\alpha = %.2f)', alpha));
    set(gca,'FontSize',11);
end

%% ======== Optional: curves at a chosen design point across alphas ========
% Pick a design point (e.g., the optimum from the mid alpha map)
ai = find(abs([results.alpha] - alpha_list(min(2,numel(alpha_list))) ) < 1e-9, 1);
if ~isempty(ai)
    % choose the (D,b) that maximizes m_max at that alpha
    [~, idx] = max(results(ai).m_max(:));
    [i_opt, j_opt]  = ind2sub(size(results(ai).m_max), idx);
    D_opt = results(ai).D_mesh(i_opt, j_opt);
    b_opt = results(ai).b_mesh(i_opt, j_opt);
    r_opt = 0.5*D_opt;

    % sweep all alphas you computed and extract Crr and torque at that (D,b)
    Crr_vs_alpha = zeros(1,numel(alpha_list));
    T_vs_alpha   = zeros(1,numel(alpha_list));
    mmax_vs_alpha= zeros(1,numel(alpha_list));
    for aidx = 1:numel(alpha_list)
        Crr_vs_alpha(aidx) = results(aidx).Crr_z(i_opt, j_opt);
        T_vs_alpha(aidx)   = results(aidx).T_wheel(i_opt, j_opt);
        mmax_vs_alpha(aidx)= results(aidx).m_max(i_opt, j_opt);
    end

    figure('Color','w');
    subplot(1,3,1);
    plot(alpha_list, Crr_vs_alpha, 'o-','LineWidth',1.5); grid on; box on;
    xlabel('\alpha = z/r'); ylabel('C_{rr}'); title(sprintf('C_{rr} vs \\alpha at D=%.2f m, b=%.2f m', D_opt, b_opt));

    subplot(1,3,2);
    plot(alpha_list, T_vs_alpha, 'o-','LineWidth',1.5); grid on; box on;
    xlabel('\alpha = z/r'); ylabel('T_{wheel} [N·m]'); title('Torque vs Sinkage');

    subplot(1,3,3);
    plot(alpha_list, mmax_vs_alpha, 'o-','LineWidth',1.5); grid on; box on;
    xlabel('\alpha = z/r'); ylabel('m_{max} [kg]'); title('m_{max} vs Sinkage');
end


