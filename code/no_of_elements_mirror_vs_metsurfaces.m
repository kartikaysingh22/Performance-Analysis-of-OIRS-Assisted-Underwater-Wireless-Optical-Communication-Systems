% ============================================================
%  Compare SNR & SE vs Number of Elements
%  - Sweep n_x from 5 to 100 (horizontal count)
%  - Keep n_m (vertical count) fixed = 13
%  - Compute both Metasurface and Mirror-Array H, then SNR & SE
%  - Plot SNR (dB) and SE on two separate figures (both curves on each)
% ============================================================

clear; clc; close all;

%% -------- Common system parameters (kept from your scripts) --------
A_PD      = 0.0001;
psi_FOV_rad = deg2rad(90);
rho       = 0.8;          % mirror reflectivity (used for mirror array)
rho_MS    = 0.8;          % metasurface reflectivity
eta       = 0.95;
T_psi     = 1;
G_psi     = 1;
c_lambda  = 0.056;
Re        = 0.5;
sigma2_i  = 1e-12;

% Lambertian
theta_half_deg = 60;
m_lambert = 2;

%% -------- Geometry & Tx/Rx positions (use the geometry you provided) --------
d_sd = 3;
x_s = 0.5;              % Tx x-position (as in your metasurface code)
x_d = x_s;              % Rx x-position aligned with Tx
z_d = d_sd/2;
y_s = 1;
y_d = 1;

S_pos = [x_s, y_s, -d_sd/2];   % Tx position
D_pos = [x_d, y_d,  z_d];      % Rx position

% patch/mirror element sizes and offsets
l_m = 0.08;   % element length (x)
w_m = 0.08;   % element width  (z)

% center offsets (use centered placement)
% We'll compute centered coordinates within loops
z_0 = 0;  % not used directly; centering uses (k-(N+1)/2)

%% -------- Sweep settings: vary n_x, keep n_m fixed = 13 --------
n_x_list = 5:5:100;     % horizontal counts (user requested 6:26)
n_m = 13;            % fixed vertical count

total_elems = n_x_list .* n_m;

% Pre-allocate results
H_MS_vals    = zeros(size(n_x_list));
H_array_vals = zeros(size(n_x_list));

fprintf('Running sweep: n_x = %d .. %d (n_m fixed = %d)\n', n_x_list(1), n_x_list(end), n_m);

% small eps to avoid singular distance (if any element exactly at Tx/Rx)
eps_dist = 1e-12;

%% -------- Loop over n_x: compute H for metasurface and mirror array --------
for idx = 1:numel(n_x_list)
    n_x = n_x_list(idx);

    %% ---- Metasurface (fixed patch normal) ----
    H_MS = 0;
    for k = 1:n_x
        % center x positions so OIRS is approximately centered at x=0
        x_kl = (k - (n_x + 1)/2) * l_m;
        for l = 1:n_m
            z_kl = (l - (n_m + 1)/2) * w_m;
            R_kl = [x_kl, 0, z_kl];

            d_sr = norm(R_kl - S_pos);
            d_rd = norm(D_pos - R_kl);
            if (d_sr < eps_dist) || (d_rd < eps_dist)
                continue;
            end

            V_SR = R_kl - S_pos;
            V_RD = D_pos - R_kl;
            V_DS = S_pos - D_pos;

            % Fixed metasurface patch normal
            n_kl = [0, 1, 0];

            % receiver FOV check: angle between V_DS and V_DR (-V_RD)
            V_DR = -V_RD;
            theta_D_R_kl = atan2(norm(cross(V_DS, V_DR)), dot(V_DS, V_DR));

            G_O_kl = 0;
            if (theta_D_R_kl <= psi_FOV_rad) && (theta_D_R_kl >= 0)
                % theta_S_R_kl: angle between -V_DS and V_SR
                theta_S_R_kl = atan2(norm(cross(-V_DS, V_SR)), dot(-V_DS, V_SR));
                c1 = (cos(theta_S_R_kl))^m_lambert;    % Lambertian Tx->patch

                c2 = cos(theta_D_R_kl);                % Rx irradiance term (patch->Rx)

                % angle between patch normal and V_RD (patch->Rx)
                theta_MS_R_kl_D = atan2(norm(cross(n_kl, V_RD)), dot(n_kl, V_RD));
                c3 = cos(theta_MS_R_kl_D);

                % angle between patch normal and -V_SR (incoming ray)
                theta_MS_R_kl_S = atan2(norm(cross(n_kl, -V_SR)), dot(n_kl, -V_SR));
                c4 = cos(theta_MS_R_kl_S);

                numerator   = (rho_MS * eta * (m_lambert + 1) * A_PD);
                denominator = 2*pi * ((d_sr + d_rd)^2);
                G_O_kl = (numerator / denominator) * c1 * c2 * c3 * c4 * T_psi * G_psi;
            end

            h_p = exp(-c_lambda * (d_sr + d_rd));   % attenuation
            H_MS = H_MS + G_O_kl * h_p;
        end
    end
    H_MS_vals(idx) = H_MS;

    %% ---- Mirror array (normals from law of reflection) ----
    H_arr = 0;
    for k = 1:n_x
        x_kl = (k - (n_x + 1)/2) * l_m;
        for l = 1:n_m
            z_kl = (l - (n_m + 1)/2) * w_m;
            R_kl = [x_kl, 0, z_kl];

            d_sr = norm(R_kl - S_pos);
            d_rd = norm(D_pos - R_kl);
            if (d_sr < eps_dist) || (d_rd < eps_dist)
                continue;
            end

            V_SR = R_kl - S_pos;
            V_RD = D_pos - R_kl;
            V_DS = S_pos - D_pos;

            u_in  = V_SR / norm(V_SR);
            u_out = V_RD / norm(V_RD);

            % Mirror element normal from law of reflection
            n_kl = ((-u_in) + u_out) / sqrt(2 + 2*((-u_in) * u_out'));

            % receiver FOV check
            V_DR = -V_RD;
            psi_kl_d = atan2(norm(cross(V_DS, V_DR)), dot(V_DS, V_DR));

            G_O_kl = 0;
            if (psi_kl_d <= psi_FOV_rad) && (psi_kl_d >= 0)
                % incidence angle at mirror between -V_SR and n_kl
                psi_s_kl = atan2(norm(cross(-V_SR, n_kl)), dot(-V_SR, n_kl));
                c2 = cos(psi_s_kl);

                % Lambertian term (Tx -> patch direction)
                phi_s_kl = atan2(norm(cross(V_SR, -V_DS)), dot(V_SR, -V_DS));
                c1 = (cos(phi_s_kl))^m_lambert;

                % irradiance term (use angle between -V_SR and normal)
                phi_kl_d = atan2(norm(cross(-V_SR, n_kl)), dot(-V_SR, n_kl));
                c3 = cos(phi_kl_d);

                % Rx incidence term
                c4 = cos(psi_kl_d);

                numerator   = (rho * eta * (m_lambert + 1) * A_PD);
                denominator = 2*pi * ((d_sr + d_rd)^2);
                G_O_kl = (numerator/denominator) * c1 * c2 * c3 * c4 * T_psi * G_psi;
            end

            h_p = exp(-c_lambda * (d_sr + d_rd));
            H_arr = H_arr + G_O_kl * h_p;
        end
    end
    H_array_vals(idx) = H_arr;

    % debug print for a particular geometry if desired
    if n_x == 13
        fprintf('DEBUG: n_x = 13 -> H_MS = %.6e, H_array = %.6e\n', H_MS, H_arr);
    end
end

%% ------- Compute SNR & SE at reference transmit power -------
Pt_dBm_ref = 25;                     % reference transmit power (dBm)
Pt_W_ref   = 1e-3 * 10^(Pt_dBm_ref/10);

P_R_MS    = Re .* Pt_W_ref .* H_MS_vals;
P_R_array = Re .* Pt_W_ref .* H_array_vals;

SNR_MS    = (P_R_MS).^2 ./ sigma2_i;
SNR_array = (P_R_array).^2 ./ sigma2_i;

SNRdB_MS    = 10 * log10(max(SNR_MS, 1e-30));
SNRdB_array = 10 * log10(max(SNR_array, 1e-30));

SE_MS    = log2(1 + SNR_MS);
SE_array = log2(1 + SNR_array);

%% ------- Plots: SNR (dB) for both on same figure -------
figure;
plot(total_elems, SNRdB_array, '-s','LineWidth',1.8,'MarkerSize',6); hold on;
plot(total_elems, SNRdB_MS,    '-o','LineWidth',1.8,'MarkerSize',6);
grid on;
xlabel('Number of elements (n_x \times n_m)');
ylabel('SNR (dB)');
title(sprintf('SNR vs Number of Elements (n_m = %d, P_t = %d dBm)', n_m, Pt_dBm_ref));
legend('Mirror Array','Metasurface','Location','best');
xlim([min(total_elems)-1 max(total_elems)+1]);

%% ------- Plots: SE for both on same figure -------
figure;
plot(total_elems, SE_array, '-s','LineWidth',1.8,'MarkerSize',6); hold on;
plot(total_elems, SE_MS,    '-o','LineWidth',1.8,'MarkerSize',6);
grid on;
xlabel('Number of elements (n_x \times n_m)');
ylabel('Spectral Efficiency (bits/s/Hz)');
title(sprintf('SE vs Number of Elements (n_m = %d, P_t = %d dBm)', n_m, Pt_dBm_ref));
legend('Mirror Array','Metasurface','Location','best');
xlim([min(total_elems)-1 max(total_elems)+1]);

%% ------- Optional: print a compact table -------
fprintf('\nSummary (n_x, total elems, H_array, H_MS, SNR_dB_array, SNR_dB_MS, SE_array, SE_MS):\n');
for i = 1:numel(n_x_list)
    fprintf('n_x=%2d, elems=%3d, H_arr=%.3e, H_MS=%.3e, SNRdB_arr=%.3f, SNRdB_MS=%.3f, SE_arr=%.4f, SE_MS=%.4f\n', ...
        n_x_list(i), total_elems(i), H_array_vals(i), H_MS_vals(i), SNRdB_array(i), SNRdB_MS(i), SE_array(i), SE_MS(i));
end
