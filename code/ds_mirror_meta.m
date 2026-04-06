% ============================================================
% SE vs Tx--Rx Distance d_sd (Mirror Array vs Metasurface)
% Sweep d_sd from 0 to 20 m, fix Pt = 25 dBm
% Both OIRS types computed and SE plotted on same figure
% ============================================================

clear; clc; close all;

%% -------------------- 1. System Parameters --------------------
A_o   = 1;            % total OIRS area (m^2)
n_elem2 = 169;        % total elements (13x13)
l_m   = 0.08;         % element length (m)
w_m   = 0.08;         % element width (m)

% optics / receiver
A_PD  = 1e-4;
psi_FOV_deg = 90;
psi_FOV_rad = deg2rad(psi_FOV_deg);
rho = 0.8;            % mirror array reflection coefficient
rho_MS = 0.8;         % metasurface reflection coefficient
eta = 0.95;
T_psi = 1;
G_psi = 1;
c_lambda = 0.056;     % attenuation (m^-1)
Re = 0.5;             % photodiode responsivity
sigma2_i = 1e-12;     % noise variance (A^2)

% Lambertian
theta_half_deg = 60;
m_lambert = 2;

%% -------------------- 2. Geometry / tiling --------------------
% fixed horizontal positions (x)
x_s = 0.5;
x_d = x_s;

% vertical depths of Tx/Rx relative to OIRS plane (y)
y_s = 1;   % Tx depth
y_d = 1;   % Rx depth (kept constant for geometry)

% OIRS tiling/pitches
x_0 = -0.56;
z_0 = -0.56;
pitch_x = l_m;
pitch_z = w_m;

% number of elements per dimension
n_m = sqrt(n_elem2);   % 13
n_p = 13;              % metasurface patches per dimension

% fixed metasurface patch normal
n_kl_fixed = [0, 1, 0];

%% -------------------- 3. Sweep settings --------------------
Pt_dBm_fixed = 30;                                    % fixed transmit power (dBm)
Pt_W_fixed   = 1e-3 * 10^(Pt_dBm_fixed/10);          % Watts

d_sd_vec = linspace(0, 20, 201);    % sweep from 0 to 20 m
eps_d = 1e-3;                        % small offset to avoid singular geometry at d_sd==0

H_array_vec = zeros(size(d_sd_vec));
H_meta_vec  = zeros(size(d_sd_vec));

SNRdB_array_vec = zeros(size(d_sd_vec));
SNRdB_meta_vec  = zeros(size(d_sd_vec));

SE_array_vec = zeros(size(d_sd_vec));
SE_meta_vec  = zeros(size(d_sd_vec));

fprintf('Starting sweep d_sd = %.2f .. %.2f m (Pt = %d dBm)\n', d_sd_vec(1), d_sd_vec(end), Pt_dBm_fixed);

%% -------------------- 4. Sweep loop --------------------
for idx = 1:numel(d_sd_vec)
    d_sd = d_sd_vec(idx);
    if d_sd == 0
        d_use = eps_d;
    else
        d_use = d_sd;
    end

    % place Tx and Rx symmetrically in z about 0 (S at -d/2, D at +d/2)
    S_pos = [x_s, y_s, -d_use/2];
    D_pos = [x_d, y_d,  d_use/2];

    %% ---- Mirror Array (normals from reflection law) ----
    H_array = 0;
    for k = 1:n_m
        x_kl = x_0 + pitch_x * k;
        for l = 1:n_m
            z_kl = z_0 + pitch_z * l;
            y_kl = 0;
            R_kl = [x_kl, y_kl, z_kl];

            d_sr = norm(R_kl - S_pos);
            d_rd = norm(D_pos - R_kl);

            % avoid degenerate cases
            if (d_sr < 1e-12) || (d_rd < 1e-12)
                continue;
            end

            V_SR = R_kl - S_pos;
            V_RD = D_pos - R_kl;
            V_DS = S_pos - D_pos;

            u_in  = V_SR / norm(V_SR);
            u_out = V_RD / norm(V_RD);

            % Mirror element normal from law of reflection (unit)
            n_kl = ((-u_in) + u_out) / sqrt(2 + 2*((-u_in) * u_out'));

            % receiver FOV check: angle between V_DS and V_DR
            V_DR = -V_RD;
            psi_kl_d = atan2(norm(cross(V_DS, V_DR)), dot(V_DS, V_DR));

            G_O_kl = 0;
            if (psi_kl_d <= psi_FOV_rad) && (psi_kl_d >= 0)
                % angle between -V_SR and n_kl (incidence)
                psi_s_kl = atan2(norm(cross(-V_SR, n_kl)), dot(-V_SR, n_kl));
                c2 = cos(psi_s_kl);

                % lambertian term: angle between V_SR and -V_DS
                phi_s_kl = atan2(norm(cross(V_SR, -V_DS)), dot(V_SR, -V_DS));
                c1 = (cos(phi_s_kl))^m_lambert;

                % irradiance geometry term (use original expression)
                e1 = [1, 0, 0]; e3 = [0, 0, 1];
                beta_kl  = asin(n_kl*e3');
                gamma_kl = asin((n_kl*e1')/max(cos(beta_kl),1e-12));
                c3 = (x_kl - x_d)*sin(beta_kl)*cos(gamma_kl)/d_rd + abs(y_kl - y_d)*cos(beta_kl)*cos(gamma_kl)/d_rd + abs(z_kl - (d_use/2))*sin(gamma_kl)/d_rd;

                % receiver term
                c4 = cos(psi_kl_d);

                numerator = rho * eta * (m_lambert + 1) * A_PD;
                denominator = 2*pi * (d_sr + d_rd)^2;
                G_O_kl = (numerator/denominator) * c1 * c2 * c3 * c4 * T_psi * G_psi;
            end

            % attenuation over S->R->D
            h_p = exp(-c_lambda * (d_sr + d_rd));

            H_array = H_array + G_O_kl * h_p;
        end
    end

    %% ---- Metasurface (fixed patch normals n_kl_fixed) ----
    H_meta = 0;
    for k = 1:n_p
        x_kl = x_0 + pitch_x * k;
        for l = 1:n_p
            z_kl = z_0 + pitch_z * l;
            y_kl = 0;
            R_kl = [x_kl, y_kl, z_kl];

            d_sr = norm(R_kl - S_pos);
            d_rd = norm(D_pos - R_kl);

            if (d_sr < 1e-12) || (d_rd < 1e-12)
                continue;
            end

            V_SR = R_kl - S_pos;
            V_RD = D_pos - R_kl;
            V_DS = S_pos - D_pos;

            n_kl = n_kl_fixed;

            V_DR = -V_RD;
            psi_d = atan2(norm(cross(V_DS, V_DR)), dot(V_DS, V_DR));

            G_O_kl = 0;
            if (psi_d <= psi_FOV_rad) && (psi_d >= 0)
                theta_SR = atan2(norm(cross(-V_DS, V_SR)), dot(-V_DS, V_SR));
                c1 = (cos(theta_SR))^m_lambert;

                c2 = cos(psi_d);

                theta_MS_R_kl_D = atan2(norm(cross(n_kl, V_RD)), dot(n_kl, V_RD));
                c3 = cos(theta_MS_R_kl_D);

                theta_MS_R_kl_S = atan2(norm(cross(n_kl, -V_SR)), dot(n_kl, -V_SR));
                c4 = cos(theta_MS_R_kl_S);

                numerator = rho_MS * eta * (m_lambert + 1) * A_PD;
                denominator = 2*pi * (d_sr + d_rd)^2;
                G_O_kl = (numerator/denominator) * c1 * c2 * c3 * c4 * T_psi * G_psi;
            end

            h_p = exp(-c_lambda * (d_sr + d_rd));
            H_meta = H_meta + G_O_kl * h_p;
        end
    end

    %% ---- Store channels ----
    H_array_vec(idx) = H_array;
    H_meta_vec(idx)  = H_meta;

    %% ---- Received power, SNR, SE ----
    P_R_array = Re * Pt_W_fixed * H_array;
    P_R_meta  = Re * Pt_W_fixed * H_meta;

    SNR_array = (P_R_array.^2) / sigma2_i;
    SNR_meta  = (P_R_meta.^2)  / sigma2_i;

    % avoid log of zero
    SNRdB_array_vec(idx) = 10*log10(max(SNR_array, 1e-30));
    SNRdB_meta_vec(idx)  = 10*log10(max(SNR_meta,  1e-30));

    % Spectral efficiency choice: use SE = log2(1 + SNR)
    % (If you prefer the other formula used earlier, replace accordingly.)
    SE_array_vec(idx) = log2(1 + exp(1) / (2 * pi) * SNR_array);
    SE_meta_vec(idx)  = log2(1 + exp(1) / (2 * pi) * SNR_meta);
end

%% -------------------- 5. Plot SE vs d_sd (both on same graph) --------------------
figure;
plot(d_sd_vec, SE_array_vec, '-s', 'LineWidth', 1.8, 'MarkerSize',5); hold on;
plot(d_sd_vec, SE_meta_vec,  '-o', 'LineWidth', 1.8, 'MarkerSize',5);
grid on;
xlabel('Tx--Rx separation d_{sd} (m)');
ylabel('Spectral Efficiency (bits/s/Hz)');
title(sprintf('SE vs Tx--Rx Distance (P_t = %d dBm)', Pt_dBm_fixed));
legend('Mirror Array','Metasurface','Location','best');
xlim([min(d_sd_vec) max(d_sd_vec)]);


%% -------------------- 7. Summary prints --------------------
[~,iAmax] = max(SE_array_vec);
[~,iMmax] = max(SE_meta_vec);
fprintf('\nMirror Array peak SE = %.4f bits/s/Hz at d_sd = %.2f m (H = %.3e)\n', SE_array_vec(iAmax), d_sd_vec(iAmax), H_array_vec(iAmax));
fprintf('Metasurface peak SE   = %.4f bits/s/Hz at d_sd = %.2f m (H = %.3e)\n\n', SE_meta_vec(iMmax), d_sd_vec(iMmax), H_meta_vec(iMmax));
