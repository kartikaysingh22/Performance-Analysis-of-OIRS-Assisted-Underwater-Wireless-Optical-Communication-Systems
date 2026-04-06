% ============================================================
% SE vs Receiver z (Mirror Array vs Metasurface) - Combined Plot
% - Fix Pt = 20 dBm
% - Tx z = -3 m (fixed)
% - Sweep Rx z from -3 to 10 m (you can change the range)
% - Both OIRS types computed and SE plotted on same figure
% ============================================================

clear; clc; close all;

%% -------------------- 1. System Parameters --------------------
A_o      = 1;            % total OIRS area (m^2) [unused directly]
n_elem2  = 169;          % total elements (13x13)
l_m      = 0.08;         % element length (m)
w_m      = 0.08;         % element width (m)
n_p      = 13;           % metasurface patches per side

% optics / receiver
A_PD     = 1e-4;
psi_FOV_deg = 90;
psi_FOV_rad = deg2rad(psi_FOV_deg);
rho      = 0.8;          % mirror array reflectivity
rho_MS   = 0.8;          % metasurface reflectivity
eta      = 0.95;
T_psi    = 1;
G_psi    = 1;
c_lambda = 0.056;        % attenuation (m^-1)
Re       = 0.5;          % photodiode responsivity (A/W)
sigma2_i = 1e-12;        % noise variance (A^2)

% Lambertian
theta_half_deg = 60;
m_lambert = 2;

%% -------------------- 2. Geometry / tiling --------------------
% fixed horizontal positions (x)
x_s = 0.5;
x_d = x_s;

% depths (y)
y_s = 1;   % Tx depth from OIRS (m)
y_d = 1;   % Rx depth from OIRS (m)

% Tx z fixed as requested
z_s = -1.5;                 % Tx z (m)
S_pos_base = [x_s, y_s, z_s];

% OIRS tiling/pitches and offsets
x_0 = -0.56;
z_0 = -0.56;
pitch_x = l_m;
pitch_z = w_m;

% mirror array elements per side
n_m = sqrt(n_elem2);

% fixed metasurface normal
n_kl_fixed = [0, 1, 0];

%% -------------------- 3. Sweep settings --------------------
Pt_dBm_fixed = 25;                           % fixed transmit power (dBm)
Pt_W_fixed   = 1e-3 * 10^(Pt_dBm_fixed/10); % W

z_d_vec = 0:0.2:10;   % sweep Rx z positions (from -3 to 10 m)
% (you said sweep Rx z from 0 to 10 in earlier snippet; I include -3->10 so Rx
% can coincide with Tx z; change the vector as you wish)

SE_meta_vec  = zeros(size(z_d_vec));
SE_array_vec = zeros(size(z_d_vec));
SNRdB_meta   = zeros(size(z_d_vec));
SNRdB_array  = zeros(size(z_d_vec));
H_meta_vec   = zeros(size(z_d_vec));
H_array_vec  = zeros(size(z_d_vec));

fprintf('Sweep: z_d = %.2f .. %.2f m, Pt = %d dBm\n', z_d_vec(1), z_d_vec(end), Pt_dBm_fixed);

%% -------------------- 4. Sweep loop --------------------
for idx = 1:numel(z_d_vec)
    z_d = z_d_vec(idx);
    S_pos = S_pos_base;
    D_pos = [x_d, y_d, z_d];

    %% ---- Metasurface  ----
    H_meta = 0;
    for k = 1:n_p
        x_kl = x_0 + pitch_x * k;
        for l = 1:n_p
            z_kl = z_0 + pitch_z * l;
            R_kl = [x_kl, 0, z_kl];

            d_sr = norm(R_kl - S_pos);
            d_rd = norm(D_pos - R_kl);
            if (d_sr < 1e-12) || (d_rd < 1e-12)
                continue;
            end

            V_SR = R_kl - S_pos;
            V_RD = D_pos - R_kl;
            V_DS = S_pos - D_pos;

            % Fixed patch normal
            n_kl = n_kl_fixed;

            V_DR = -V_RD;
            psi_d = atan2(norm(cross(V_DS, V_DR)), dot(V_DS, V_DR)); % angle for FOV check

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

    %% ---- Mirror Array ----
    H_array = 0;
    for k = 1:n_m
        x_kl = x_0 + pitch_x * k;
        for l = 1:n_m
            z_kl = z_0 + pitch_z * l;
            R_kl = [x_kl, 0, z_kl];

            d_sr = norm(R_kl - S_pos);
            d_rd = norm(D_pos - R_kl);
            if (d_sr < 1e-12) || (d_rd < 1e-12)
                continue;
            end

            V_SR = R_kl - S_pos;
            V_RD = D_pos - R_kl;
            V_DS = S_pos - D_pos;

            u_in  = V_SR / norm(V_SR);
            u_out = V_RD / norm(V_RD);

            % mirror normal by reflection law
            n_kl = ((-u_in) + u_out) / sqrt(2 + 2*((-u_in) * u_out'));

            V_DR = -V_RD;
            psi_kl_d = atan2(norm(cross(V_DS, V_DR)), dot(V_DS, V_DR));

            G_O_kl = 0;
            if (psi_kl_d <= psi_FOV_rad) && (psi_kl_d >= 0)
                psi_s_kl = atan2(norm(cross(-V_SR, n_kl)), dot(-V_SR, n_kl));
                c2 = cos(psi_s_kl);

                phi_s_kl = atan2(norm(cross(V_SR, -V_DS)), dot(V_SR, -V_DS));
                c1 = (cos(phi_s_kl))^m_lambert;

                % use angle between -V_SR and n_kl for irradiance term (safer numerics)
                phi_kl_d = atan2(norm(cross(-V_SR, n_kl)), dot(-V_SR, n_kl));
                c3 = cos(phi_kl_d);

                c4 = cos(psi_kl_d);

                numerator = rho * eta * (m_lambert + 1) * A_PD;
                denominator = 2*pi * (d_sr + d_rd)^2;
                G_O_kl = (numerator/denominator) * c1 * c2 * c3 * c4 * T_psi * G_psi;
            end

            h_p = exp(-c_lambda * (d_sr + d_rd));
            H_array = H_array + G_O_kl * h_p;
        end
    end

    %% ---- store H and compute SE ----
    H_meta_vec(idx)  = H_meta;
    H_array_vec(idx) = H_array;

    % Received electrical current (A)
    P_R_meta  = Re * Pt_W_fixed * H_meta;
    P_R_array = Re * Pt_W_fixed * H_array;

    SNR_meta = (P_R_meta.^2) / sigma2_i;
    SNR_arr  = (P_R_array.^2) / sigma2_i;

    SNRdB_meta(idx)  = 10*log10(max(SNR_meta, 1e-30));
    SNRdB_array(idx) = 10*log10(max(SNR_arr,  1e-30));

    % Spectral efficiency (consistent formula for both)
    SE_meta_vec(idx)  = log2(1 + exp(1) / (2 * pi) * SNR_meta);
    SE_array_vec(idx) = log2(1 + exp(1) / (2 * pi) * SNR_arr);
end

%% -------------------- 5. Plot SE (both on same figure) --------------------
figure;
plot(z_d_vec, SE_array_vec, '-s', 'LineWidth', 1.8, 'MarkerSize',5); hold on;
plot(z_d_vec, SE_meta_vec,  '-o', 'LineWidth', 1.8, 'MarkerSize',5);
grid on;
xlabel('Receiver z position (m)');
ylabel('Spectral Efficiency (bits/s/Hz)');
title(sprintf('SE vs Receiver z (P_t = %d dBm, Tx z = %.1f m)', Pt_dBm_fixed, z_s));
legend('Mirror Array','Metasurface','Location','best');
xlim([min(z_d_vec) max(z_d_vec)]);

%% -------------------- 6. additional diagnostics --------------------
% Plot SNRs for reference
figure;
plot(z_d_vec, SNRdB_array, '-s', 'LineWidth', 1.6); hold on;
plot(z_d_vec, SNRdB_meta,  '-o', 'LineWidth', 1.6);
grid on;
xlabel('Receiver z position (m)');
ylabel('SNR (dB)');
title('SNR vs Receiver z');
legend('Mirror Array','Metasurface','Location','best');

% Quick printout of peaks
[~, iA] = max(SE_array_vec);
[~, iM] = max(SE_meta_vec);
fprintf('\nMirror Array peak SE = %.4f bits/s/Hz at z = %.2f m (H = %.3e)\n', SE_array_vec(iA), z_d_vec(iA), H_array_vec(iA));
fprintf('Metasurface peak SE   = %.4f bits/s/Hz at z = %.2f m (H = %.3e)\n\n', SE_meta_vec(iM), z_d_vec(iM), H_meta_vec(iM));
