%% Combined: Mirror Array vs Metasurface vs Plane Mirror (OIRS)
clear; clc; close all;

%% -------------------- Common System Parameters --------------------
A_o   = 1;              % Total reflective surface area (m^2)
n_elem2 = 169;          % Mirror array: total number of mirror elements (13×13)
l_m   = 0.08;           % Mirror element length (m)
w_m   = 0.08;           % Mirror element width (m)
y_s   = 1;              % Tx depth relative to OIRS (m)
y_d   = 1;              % Rx depth relative to OIRS (m)
n_m = 13;               % mirror elements per dimension (13)

% Metasurface parameters
n_p = 13;               % Number of metasurface patches per dimension (13×13)

% EGG turbulence / physical constants (not used here)
omega = 0.1665;
lambda_exp = 0.1207;
a_gg = 0.1559;
b_gg = 1.5216;
c_gg = 22.8754;

% Receiver / optics
A_PD  = 0.0001;         % Photodiode area (m^2)
psi_FOV_deg = 90;      
psi_FOV_rad = deg2rad(psi_FOV_deg);
mu = 1.33;             
rho = 0.8;              % Reflection coefficient — mirror array
rho_MS = 0.8;           % Reflection coefficient — metasurface
eta = 0.95;             
T_psi = 1;              
G_psi = 1;              
c_lambda = 0.056;       % Attenuation coefficient (m^-1)
Re = 0.5;               % Photodiode responsivity
sigma2_i = 1e-12;       % Noise variance (A^2)
c_water = 2.25e8;       

% Lambertian source
theta_half_deg = 60;  
m_lambert = 2;         

%% -------------------- Geometry --------------------
d_sd = 3;                       % Tx–Rx separation (m)
x_s = 0.5;                      
x_d = x_s;                      
z_d = d_sd/2;                   

S_pos = [x_s, y_s, -d_sd/2];    % Tx coordinates
D_pos = [x_s, y_d, z_d];        % Rx coordinates

%% -------------------- 1) Compute H for Mirror Array OIRS --------------------
x_0 = -0.56;
z_0 = -0.56;
H_array = 0;

for k = 1:n_m
    x_kl = x_0 + l_m * k;
    for l = 1:n_m
        z_kl = z_0 + w_m * l;
        y_kl = 0;
        R_kl = [x_kl, y_kl, z_kl];

        d_sr = norm(R_kl - S_pos);
        d_rd = norm(D_pos - R_kl);

        V_SR = R_kl - S_pos;
        V_RD = D_pos - R_kl;
        V_DS = S_pos - D_pos;

        u_in  = V_SR / norm(V_SR);
        u_out = V_RD / norm(V_RD);

        % Compute surface normal based on reflection law
        n_kl = ((-u_in)+u_out)/sqrt(2+2*((-u_in) * u_out'));

        e1 = [1, 0, 0];
        e3 = [0, 0, 1];

        beta_kl = asin(n_kl*e3');
        gamma_kl = asin((n_kl*e1')/cos(beta_kl));

        V_DR = -V_RD;
        psi_kl_d = atan2(norm(cross(V_DS,V_DR)),dot(V_DS,V_DR));

        G_O_kl = 0;

        if (psi_kl_d <= psi_FOV_rad) && (psi_kl_d >= 0)
            psi_s_kl = atan2(norm(cross(-V_SR,n_kl)),dot(-V_SR,n_kl));
            c2 = cos(psi_s_kl);

            phi_s_kl = atan2(norm(cross(V_SR,-V_DS)),dot(V_SR,-V_DS));
            c1 = (cos(phi_s_kl))^m_lambert;

            c3 = (x_kl - x_d)*sin(beta_kl)*cos(gamma_kl)/d_rd + abs(y_kl-y_d)*cos(beta_kl)*cos(gamma_kl)/d_rd + abs(z_kl - z_d)*sin(gamma_kl)/d_rd;
            c4 = cos(psi_kl_d);

            numerator = (rho* eta * (m_lambert+1)*A_PD);
            denominator = 2*pi*((d_sr+d_rd)^2);
            G_O_kl = (numerator/denominator)*c1*c2*c3*c4*T_psi*G_psi;
        end

        h_p = exp(-c_lambda * (d_sr + d_rd));
        H_array = H_array + G_O_kl * h_p;
    end
end
fprintf('H_array = %.6e\n', H_array);

%% -------------------- 2) Compute H for Metasurface OIRS --------------------
H_meta = 0;
x_0 = -0.56;
z_0 = -0.56;
for k = 1:n_p
    x_kl = x_0 + l_m * k;
    for l = 1:n_p
        z_kl = z_0 + 0.08 * l;
        y_kl = 0;
        R_kl = [x_kl, y_kl, z_kl];

        d_sr = norm(R_kl - S_pos);
        d_rd = norm(D_pos - R_kl);

        V_SR = R_kl - S_pos;
        V_RD = D_pos - R_kl;
        V_DS = S_pos - D_pos;
        u_in  = V_SR/norm(V_SR);
        u_out = V_RD / norm(V_RD);
        n_kl = [0,1,0];

        V_DR = -V_RD;
        theta_D_R_kl = atan2(norm(cross(V_DS,V_DR)),dot(V_DS,V_DR));

        G_O_kl = 0;
        if ((theta_D_R_kl<=psi_FOV_rad) && (theta_D_R_kl>=0))
            theta_S_R_kl = atan2(norm(cross(-V_DS,V_SR)),dot(-V_DS,V_SR));
            c1 = cos(theta_S_R_kl)^m_lambert;

            c2 = cos(theta_D_R_kl);

            theta_MS_R_kl_D = atan2(norm(cross(n_kl,V_RD)),dot(n_kl,V_RD));
            c3 = cos(theta_MS_R_kl_D);

            theta_MS_R_kl_S = atan2(norm(cross(n_kl,-V_SR)),dot(n_kl,-V_SR));
            c4 = cos(theta_MS_R_kl_S);

            numerator = (rho_MS* eta* (m_lambert+1)*A_PD);
            denominator = 2*pi*((d_sr+d_rd)^2) * norm(n_kl) * norm(u_out);
            G_O_kl = (numerator/denominator)*c1*c2*c3*c4*T_psi*G_psi;
        end

        h_p = exp(-c_lambda * (d_sr + d_rd));
        H_meta = H_meta + G_O_kl * h_p;
    end
end
fprintf('H_meta = %.6e\n', H_meta);

%% -------------------- 3) Compute H for Plane Mirror (OIRS) --------------------
% Using the provided plane OIRS block; note some geometry params differ in that block.
% I'll keep most parameters consistent but use the plane block geometry (y_s_plane,y_d_plane,d_sd_plane)
% so plane mirror H is computed per the provided code.

% Plane-specific geometry (from provided plane mirror code)
y_s_plane = 2; y_d_plane = 2; d_sd_plane = 6;
S_pos_plane = [1, y_s_plane, -d_sd_plane/2];
D_pos_plane = [1, y_d_plane, d_sd_plane/2];

n_m_plane = sqrt(n_elem2);  % should be 13
H_plane = 0;
x_0 = -0.56;
z_0 = -0.56;
for k = 1:n_m_plane
    x_kl = x_0 + 0.08 * k;
    for l = 1:n_m_plane
        z_kl = z_0 + 0.08 * l;
        y_kl = 0;
        R_kl = [x_kl, 0, z_kl];

        d_sr = norm(R_kl - S_pos_plane);
        d_rd = norm(D_pos_plane - R_kl);

        V_SR = R_kl - S_pos_plane;
        V_RD = D_pos_plane - R_kl;
        V_DS = S_pos_plane - D_pos_plane;
        u_in  = V_SR/norm(V_SR);
        u_out = V_RD/norm(V_RD);
        n_kl = [0,1,0];

        V_DR = -V_RD;
        psi_kl_d = atan2(norm(cross(V_DS,V_DR)),dot(V_DS,V_DR));

        G_O_kl = 0;
        if ((psi_kl_d<=psi_FOV_rad) && (psi_kl_d>=0))
            psi_s_kl = atan2(norm(cross(-V_SR,n_kl)),dot(-V_SR,n_kl));
            c2 = cos(psi_s_kl);

            phi_s_kl = atan2(norm(cross(V_SR,-V_DS)),dot(V_SR,-V_DS));
            c1 = (cos(phi_s_kl))^m_lambert;

            phi_kl_d = atan2(norm(cross(-V_SR,n_kl)),dot(-V_SR,n_kl));
            c3 = cos(phi_kl_d);

            c4 = cos(psi_kl_d);

            numerator = (rho* eta * (m_lambert+1)*A_PD);
            denominator = 2*pi*((d_sr+d_rd)^2);
            G_O_kl = (numerator/denominator)*c1*c2*c3*c4*T_psi*G_psi;
        else
            G_O_kl = 0;
        end

        h_p = exp(-c_lambda*(d_sr+d_rd));
        H_plane = H_plane + G_O_kl * h_p;
    end
end
fprintf('H_plane (OIRS) = %.6e\n', H_plane);

%% -------------------- 4) Power Sweep, SNR and Spectral Efficiency --------------------
Pt_dBm = 0:3:60;                         % Power sweep (dBm)
Pt_W   = 1e-3 * 10.^(Pt_dBm/10);        % Convert to Watts

% Received power
P_R_array = Re .* Pt_W .* H_array;       
P_R_meta  = Re .* Pt_W .* H_meta;        
P_R_plane = Re .* Pt_W .* H_plane;

% Electrical SNR (as in your original script)
SNR_array = (P_R_array).^2 ./ sigma2_i;
SNR_meta  = (P_R_meta).^2  ./ sigma2_i;
SNR_plane = (P_R_plane).^2 ./ sigma2_i;

SNRdB_array = 10*log10(SNR_array );
SNRdB_meta  = 10*log10(SNR_meta  );
SNRdB_plane = 10*log10(SNR_plane );

% Use the same SE formula you used earlier for array/meta so comparison is fair:
SE_array = log2(1 + exp(1) / (2 * pi) * SNR_array);
SE_meta  = log2(1 + exp(1) / (2 * pi) * SNR_meta);
SE_plane = log2(1 + exp(1) / (2 * pi) * SNR_plane);   % note: originally plane code used log2(1+SNR)

%% -------------------- 5) Plot SNR and Spectral Efficiency (all three) --------------------
figure('Units','normalized','Position',[0.05 0.05 0.9 0.5]);

% --- SNR Plot ---
subplot(1,2,1);
plot(Pt_dBm, SNRdB_array, '-o','LineWidth',1.8); hold on;
plot(Pt_dBm, SNRdB_meta,  's-','LineWidth',1.8);
plot(Pt_dBm, SNRdB_plane, 'd-','LineWidth',1.8);
grid on;
xlabel('Transmit Power (dBm)');
ylabel('SNR (dB)');
title('SNR vs Transmit Power');
legend('Mirror Array','Metasurface','Plane Mirror (OIRS)','Location','SouthEast');

% --- Spectral Efficiency Plot ---
subplot(1,2,2);
plot(Pt_dBm, SE_array, '-o','LineWidth',1.8); hold on;
plot(Pt_dBm, SE_meta,  's-','LineWidth',1.8);
plot(Pt_dBm, SE_plane, 'd-','LineWidth',1.8);
grid on;
xlabel('Transmit Power (dBm)');
ylabel('Spectral Efficiency (bits/s/Hz)');
title('Spectral Efficiency Comparison');
legend('Mirror Array','Metasurface','Plane Mirror (OIRS)','Location','NorthWest');

%% (Optional) quick numeric output at a few transmit powers
fprintf('\nSample values at Pt = [0 30 60] dBm:\n');
idxs = [1, find(Pt_dBm==30), find(Pt_dBm==60)];
for i = idxs
    fprintf('Pt=%2d dBm -> SE_array=%.3f, SE_meta=%.3f, SE_plane=%.3f\n', ...
        Pt_dBm(i), SE_array(i), SE_meta(i), SE_plane(i));
end
