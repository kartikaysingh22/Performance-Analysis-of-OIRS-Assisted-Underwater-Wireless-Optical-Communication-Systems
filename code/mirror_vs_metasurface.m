%% =======================================================================
%   Combined Comparison: Mirror Array vs Metasurface
%   Computes and compares channel gain, SNR, and SE for:
%       1) Conventional Mirror Array OIRS
%       2) Metasurface-based OIRS
% =======================================================================

clear; clc; close all;

%% -------------------- 1. System Parameters --------------------
A_o   = 1;              % Total reflective surface area (m^2)
n_elem2 = 169;          % Mirror array: total number of mirror elements (13×13)
l_m   = 0.08;           % Mirror element length (m)
w_m   = 0.08;           % Mirror element width (m)
y_s   = 1;              % Tx depth relative to OIRS (m)
y_d   = 1;              % Rx depth relative to OIRS (m)
n_m = 13;            % mirror elements per dimension (13)

% --- Metasurface parameters ---
n_p = 13;               % Number of metasurface patches per dimension (13×13)

% --- EGG turbulence / physical constants (not used yet) ---
omega = 0.1665;
lambda_exp = 0.1207;
a_gg = 0.1559;
b_gg = 1.5216;
c_gg = 22.8754;

% --- Receiver / optics parameters ---
A_PD  = 0.0001;         % Photodiode area (m^2)
psi_FOV_deg = 90;      % Field of view (deg)
psi_FOV_rad = deg2rad(psi_FOV_deg);   % Convert FoV to radians
mu = 1.33;             % Refractive index (water)
rho = 0.8;             % Reflection coefficient — mirror array
rho_MS = 0.8;          % Reflection coefficient — metasurface
eta = 0.95;            % LED efficiency
T_psi = 1;             % Optical filter gain
G_psi = 1;             % Optical concentrator gain
c_lambda = 0.056;      % Attenuation coefficient (m^-1)
Re = 0.5;              % Photodiode responsivity
sigma2_i = 1e-12;      % Noise variance (A^2)
c_water = 2.25e8;      % Speed of light in water (m/s)

% --- Lambertian source parameters ---
theta_half_deg = 60;  
m_lambert = 2;         % Lambertian order (approx)

%% -------------------- 2. Geometry --------------------
d_sd = 3;                       % Tx–Rx separation (m)
x_s = 0.5;                        % Tx x-position
x_d = x_s;                      % Rx x-position
z_d = d_sd/2;                   % Rx z-position

S_pos = [x_s, y_s, -d_sd/2];    % Tx coordinates
D_pos = [x_s, y_d, z_d];        % Rx coordinates


%% -------------------- 3. Compute H for Mirror Array OIRS --------------------

% Starting coordinates for mirror tile grid
x_0 = -0.56;    % k will update x
z_0 = -0.56;    % l will update z     
H_array = 0;                    % Accumulated channel gain

for k = 1:n_m
    x_kl = x_0 + l_m * k;       % x-position for element (k,l)
    for l = 1:n_m
        z_kl = z_0 + w_m * l;   % z-position for element (k,l)
        y_kl = 0;               % OIRS plane at y = 0
        R_kl = [x_kl, y_kl, z_kl];  % Reflective element coordinates

        % Distances Tx→Reflector and Reflector→Rx
        d_sr = norm(R_kl - S_pos);
        d_rd = norm(D_pos - R_kl);

        % Key geometry vectors
        V_SR = R_kl - S_pos;    % Source(S) to Reflector
        V_RD = D_pos - R_kl;    % Reflector to Destination (D)
        V_DS = S_pos - D_pos;   % D to S

        % Normalize ray directions
        u_in  = V_SR / (norm(V_SR));
        u_out = V_RD / (norm(V_RD));

        % Compute surface normal based on reflection law
         n_kl = ((-u_in)+u_out)/sqrt(2+2*((-u_in) * u_out')); %%% normal vector

        e1 = [1, 0, 0];
        e3 = [0, 0, 1]; %3rd column of identity matrix


        beta_kl = asin(n_kl*e3'); 
        gamma_kl =asin((n_kl*e1')/cos(beta_kl));
        

        V_DR = -V_RD;
        psi_kl_d = atan2(norm(cross(V_DS,V_DR)),dot(V_DS,V_DR)); %% angle between V_DS and V_DR

        G_O_kl = 0;             % Initialize per-element gain

        if (psi_kl_d <= psi_FOV_rad) && (psi_kl_d >= 0)
           %%% Calculation of G_O_kl

            %calculating angle between V_RS and n_kl, ie, 
            % angle of incidence
            % psi_s_kl 
            %needed for c2
            psi_s_kl = atan2(norm(cross(-V_SR,n_kl)),dot(-V_SR,n_kl));
            c2 = cos(psi_s_kl);

            %calculating phi_s_kl
            %it is the angle between V_SR and V_SD(-V_DS)
            %needed for c1
            phi_s_kl = atan2(norm(cross(V_SR,-V_DS)),dot(V_SR,-V_DS));
            c1 = (cos(phi_s_kl))^m_lambert;

            %%calculating c3 = cos (phi_kl_d) whchi is angle of irradiance
            c3 = (x_kl - x_d)*sin(beta_kl)*cos(gamma_kl)/d_rd + abs(y_kl-y_d)*cos(beta_kl)*cos(gamma_kl)/d_rd + abs(z_kl - z_d)*sin(gamma_kl)/d_rd;
            %phi_kl_d = atan2(norm(cross(-V_SR,n_kl)),dot(-V_SR,n_kl));
            %c3 = cos(phi_kl_d);

            %c4 is cos(psi_kl_d)
            %angle between V_DR and V_DS
            c4 = cos(psi_kl_d);

            %calculating G_O_kl
            numerator = (rho* eta * (m_lambert+1)*A_PD);
            denominator = 2*pi*((d_sr+d_rd)^2);
            G_O_kl = (numerator/denominator)*c1*c2*c3*c4*T_psi*G_psi;
        end

        % Water attenuation
        h_p = exp(-c_lambda * (d_sr + d_rd));

        % Add this element's contribution
        H_array = H_array + G_O_kl * h_p;
    end
end
fprintf('H_array = %.6e\n', H_array);

%% -------------------- 4. Compute H for Metasurface OIRS --------------------
H_meta = 0;
% Starting coordinates for mirror tile grid
x_0 = -0.56;    % k will update x
z_0 = -0.56;    % l will update z 
for k = 1:n_p
    x_kl = x_0 + l_m * k;
    for l = 1:n_p
        z_kl = z_0 + 0.08 * l;
        y_kl = 0;
        R_kl = [x_kl, y_kl, z_kl];

        d_sr = norm(R_kl - S_pos);
        d_rd = norm(D_pos - R_kl);

        V_SR = R_kl - S_pos; %vector from S-pos to Reflector
        V_RD = D_pos - R_kl; %vector from D-pos to Reflector
        V_DS = S_pos - D_pos; % vector from D-pos to S-pos
        u_in  = V_SR/norm(V_SR); %unit incidence ray
        u_out =  V_RD / norm(V_RD);%unit reflected ray
        n_kl = [0,1,0]; % normal of a patch remains same for each patch
        
        V_DR = -V_RD;
        theta_D_R_kl = atan2(norm(cross(V_DS,V_DR)),dot(V_DS,V_DR)); %% angle between V_DS and V_DR

        G_O_kl = 0; %%%intilization of OIRS GAIN
        if ((theta_D_R_kl<=psi_FOV_rad) && (theta_D_R_kl>=0))
            %%% Calculation of G_O_kl

            %%% theta_S_R_kl
            %%%angle between SD and SR_kl
            theta_S_R_kl = atan2(norm(cross(-V_DS,V_SR)),dot(-V_DS,V_SR));
            c1 = cos(theta_S_R_kl)^m_lambert;

            % theta_D_R_kl already calculated
            c2 = cos(theta_D_R_kl);

            % theta_MS_R_kl_D
            % angle between MS normal n_kl and V_RD
            theta_MS_R_kl_D = atan2(norm(cross(n_kl,V_RD)),dot(n_kl,V_RD));
            c3 = cos(theta_MS_R_kl_D);

            %theta_MS_R_kl_S
            %angle bw MS normal and minus V_S_R
            theta_MS_R_kl_S = atan2(norm(cross(n_kl,-V_SR)),dot(n_kl,-V_SR));
            c4 = cos(theta_MS_R_kl_S);
 

            %calculating G_O_kl
            numerator = (rho_MS* eta* (m_lambert+1)*A_PD);
            denominator = 2*pi*((d_sr+d_rd)^2) * norm(n_kl) * norm(u_out);
            G_O_kl = (numerator/denominator)*c1*c2*c3*c4*T_psi*G_psi;
        end

        h_p = exp(-c_lambda * (d_sr + d_rd)); %%% attentuation due to path loss
        H_meta = H_meta + G_O_kl * h_p;
    end
end
fprintf('H_meta = %.6e\n', H_meta);

%% -------------------- 5. Power Sweep, SNR and Spectral Efficiency --------------------
Pt_dBm = 0:3:60;                         % Power sweep (dBm)
Pt_W   = 1e-3 * 10.^(Pt_dBm/10);         % Convert to Watts

P_R_array = Re .* Pt_W .* H_array;       % Received power (mirror array)
P_R_meta  = Re .* Pt_W .* H_meta;        % Received power (metasurface)

SNR_array = (P_R_array).^2 ./ sigma2_i;  % Electrical SNR
SNR_meta  = (P_R_meta).^2  ./ sigma2_i;

SNRdB_array = 10*log10(SNR_array );
SNRdB_meta  = 10*log10(SNR_meta  );

SE_array = log2(1 + exp(1) / (2 * pi) * SNR_array);          % Shannon capacity
SE_meta  = log2(1 + exp(1) / (2 * pi) * SNR_meta);

%% -------------------- 6. Plot SNR and Spectral Efficiency --------------------
figure('Units','normalized','Position',[0.05 0.05 0.9 0.5]);

% --- SNR Plot ---
subplot(1,2,1);
plot(Pt_dBm, SNRdB_array, '-o','LineWidth',1.8); hold on;
plot(Pt_dBm, SNRdB_meta,  's-','LineWidth',1.8);
grid on;
xlabel('Transmit Power (dBm)');
ylabel('SNR (dB)');
title('SNR vs Transmit Power');
legend('Mirror Array','Metasurface','Location','SouthEast');

% --- Spectral Efficiency Plot ---
subplot(1,2,2);
plot(Pt_dBm, SE_array, '-o','LineWidth',1.8); hold on;
plot(Pt_dBm, SE_meta,  's-','LineWidth',1.8);
grid on;
xlabel('Transmit Power (dBm)');
ylabel('Spectral Efficiency (bits/s/Hz)');
title('Spectral Efficiency Comparison');
legend('Mirror Array','Metasurface','Location','NorthWest');
