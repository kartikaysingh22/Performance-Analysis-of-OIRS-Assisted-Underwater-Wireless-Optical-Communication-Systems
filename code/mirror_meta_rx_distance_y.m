% ============================================================
% SE vs Receiver Horizontal Distance (Mirror Array vs Metasurface)
% Fixed z_d = d_sd/2, Sweep y_d from 0 to 20 m
% Pt = 30 dBm for both
% ============================================================

clear; clc; close all;

%% ===================== COMMON PARAMETERS =====================
A_PD  = 1e-4;
psi_FOV_deg = 90;
psi_FOV_rad = deg2rad(psi_FOV_deg);
rho = 0.8;
rho_MS = 0.8;
eta = 0.95;
T_psi = 1;
G_psi = 1;
c_lambda = 0.056;
Re = 0.5;
sigma2_i = 1e-12;
m_lambert = 2;

%% ===================== GEOMETRY SETTINGS =====================
d_sd = 4;
x_s = 0;
x_d = 0;
y_s = 1;
z_s = -d_sd/2;
z_d = d_sd/2;

S_pos = [x_s, y_s, z_s];

y_d_vec = 0:0.2:10;

Pt_dBm_fixed = 30;
Pt_W_fixed = 1e-3 * 10^(Pt_dBm_fixed/10);

%% ===================== MIRROR ARRAY SETTINGS =====================
n_elem2 = 169;
n_m = sqrt(n_elem2);
l_m = 0.08;
w_m = 0.08;

x_0 = -0.56;
z_0 = -0.56;

SE_mirror = zeros(size(y_d_vec));

%% ===================== MIRROR ARRAY SIMULATION =====================
for iy = 1:numel(y_d_vec)

    y_d = y_d_vec(iy);
    D_pos = [x_d, y_d, z_d];
    H_OIRS = 0;

    for k = 1:n_m
        x_kl = x_0 + l_m*k;
        for l = 1:n_m
            z_kl = z_0 + w_m*l;
            R_kl = [x_kl, 0, z_kl];

            d_sr = norm(R_kl - S_pos);
            d_rd = norm(D_pos - R_kl);

            V_SR = R_kl - S_pos;
            V_RD = D_pos - R_kl;
            V_DS = S_pos - D_pos;

            u_in  = V_SR/norm(V_SR);
            u_out = V_RD/norm(V_RD);

            n_kl = ((-u_in)+u_out)/sqrt(2+2*((-u_in)*u_out'));

            e1 = [1 0 0]; e3 = [0 0 1];
            beta  = asin(n_kl*e3');
            gamma = asin((n_kl*e1')/cos(beta));

            V_DR = -V_RD;
            psi_d = atan2(norm(cross(V_DS,V_DR)),dot(V_DS,V_DR));

            G = 0;
            if psi_d <= psi_FOV_rad
                psi_s = atan2(norm(cross(-V_SR,n_kl)),dot(-V_SR,n_kl));
                c2 = cos(psi_s);

                phi_s = atan2(norm(cross(V_SR,-V_DS)),dot(V_SR,-V_DS));
                c1 = (cos(phi_s))^m_lambert;

                c3 = (x_kl-x_d)*sin(beta)*cos(gamma)/d_rd + ...
                     abs(y_d)*cos(beta)*cos(gamma)/d_rd + ...
                     abs(z_kl-z_d)*sin(gamma)/d_rd;

                c4 = cos(psi_d);

                numerator = rho*eta*(m_lambert+1)*A_PD;
                denominator = 2*pi*(d_sr+d_rd)^2;
                G = (numerator/denominator)*c1*c2*c3*c4*T_psi*G_psi;
            end

            h_p = exp(-c_lambda*(d_sr+d_rd));
            H_OIRS = H_OIRS + G*h_p;
        end
    end

    P_R = Re * Pt_W_fixed * H_OIRS;
    SNR = (P_R^2)/sigma2_i;
    SE_mirror(iy) = log2(1 + SNR);
end

%% ===================== METASURFACE SETTINGS =====================
n_p = 13;
pitch_x = l_m;
pitch_z = w_m;

SE_meta = zeros(size(y_d_vec));

%% ===================== METASURFACE SIMULATION =====================
for iy = 1:numel(y_d_vec)

    y_d = y_d_vec(iy);
    D_pos = [x_d, y_d, z_d];
    H_MS = 0;

    for k = 1:n_p
        x_kl = x_0 + pitch_x*k;
        for l = 1:n_p
            z_kl = z_0 + pitch_z*l;
            R_kl = [x_kl, 0, z_kl];

            d_sr = norm(R_kl - S_pos);
            d_rd = norm(D_pos - R_kl);

            V_SR = R_kl - S_pos;
            V_RD = D_pos - R_kl;
            V_DS = S_pos - D_pos;

            n_kl = [0 1 0];
            V_DR = -V_RD;
            psi_d = atan2(norm(cross(V_DS,V_DR)),dot(V_DS,V_DR));

            G = 0;
            if psi_d <= psi_FOV_rad

                theta_SR = atan2(norm(cross(-V_DS,V_SR)),dot(-V_DS,V_SR));
                c1 = (cos(theta_SR))^m_lambert;

                c2 = cos(psi_d);

                c3 = cos(atan2(norm(cross(n_kl,V_RD)),dot(n_kl,V_RD)));
                c4 = cos(atan2(norm(cross(n_kl,-V_SR)),dot(n_kl,-V_SR)));

                numerator = rho_MS*eta*(m_lambert+1)*A_PD;
                denominator = 2*pi*(d_sr+d_rd)^2;
                G = (numerator/denominator)*c1*c2*c3*c4*T_psi*G_psi;
            end

            h_p = exp(-c_lambda*(d_sr+d_rd));
            H_MS = H_MS + G*h_p;
        end
    end

    P_R = Re * Pt_W_fixed * H_MS;
    SNR = (P_R^2)/sigma2_i;
    SE_meta(iy) = log2(1 + SNR);
end

%% ===================== COMBINED COMPARISON PLOT =====================
figure;
plot(y_d_vec, SE_mirror,'-s','LineWidth',1.8,'MarkerFaceColor','auto'); hold on;
plot(y_d_vec, SE_meta,  '-o','LineWidth',1.8,'MarkerFaceColor','auto');
grid on;
xlabel('Receiver vertical distance y_d (m)');
ylabel('Spectral Efficiency (bits/s/Hz)');
title('Spectral Efficiency vs Receiver Vertical Distance (Pt = 30 dBm)');
legend('Mirror Array','Metasurface','Location','best');

