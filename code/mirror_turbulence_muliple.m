% ============================================================
%   OIRS-assisted UWOC Simulation with EGG Turbulence
%   MIRROR ARRAY (Rx NORMAL IS HORIZONTAL)
% ============================================================
clear; clc; close all;

%% -------------------- 1. System Parameters --------------------
A_o   = 1;          % Total OIRS area (m^2)
n_elem2 = 169;      % Total number of mirror elements (n_m^2)
l_m   = 0.08;       % Mirror element length (m)
w_m   = 0.08;       % Mirror element width (m)
y_s   = 1;         % Tx depth from OIRS (m)
y_d   = 1;         % Rx depth from OIRS (m)

% --- EGG turbulence parameters (values taken from table) ---
omega = 0.1665;        % mixture weight
lambda_exp = 0.1207;   % exponential parameter (scale)
a_gg = 0.1559;         % gamma shape a
b_gg = 1.5216;         % gen-gamma b (scale)
c_gg = 22.8754;        % gen-gamma c (power)

% --- Physical / receiver constants ---
A_PD  = 0.0001;       % photodiode area (m^2)
psi_FOV_deg = 90;     % receiver FoV (deg)
psi_FOV_rad = deg2rad(psi_FOV_deg); % conversion into radians
mu = 1.33;            % refractive index (water)
rho = 0.8;            % mirror reflection coefficient
eta = 0.95;           % LED efficiency
T_psi = 1;            % optical filter gain
G_psi = 1;
c_lambda = 0.056;     % attenuation coefficient (m^-1)
Re = 0.5;             % responsivity of photodiode (A/W)
sigma2_i = 1e-12;     % total noise variance (A^2)
c_water = 2.25e8;     % light speed in water (m/s)

% --- LED Lambertian source parameters ---
theta_half_deg = 60;  
m_lambert = 2;        %m_lambert= -log(2)/log(cosd(theta_half_deg)) fomula
%%%% lamberterian order


%% -------------------- 2. Geometry --------------------
n_m = sqrt(n_elem2);
d_sd = 4;
% Tx-Rx separation (m)
x_s = 0.5; 
x_d = x_s;
z_d = d_sd/2;
S_pos = [x_s, y_s, -d_sd/2];    % Tx position
D_pos = [x_s, y_d, z_d];       % Rx position


%% -------------------- 3. Compute OIRS Channel Gain (deterministic H_OIRS) ----
fprintf('Computing per-element OIRS gains...\n');
H_OIRS = 0;
count_nonzero = 0;
total_length = n_m * l_m;
total_width  = n_m * w_m;
x_0 = -0.56; %% k will update x
z_0 = -0.56; %% l will update z


for k = 1:n_m
    x_kl = x_0 + 0.08 * k;
    for l = 1:n_m
        z_kl = z_0 + 0.08 * l;
        y_kl = 0;
        R_kl = [x_kl, 0, z_kl];

        d_sr = norm(R_kl - S_pos);
        d_rd = norm(D_pos - R_kl);

        V_SR = R_kl - S_pos; %vector from S-pos to Reflector
        V_RD = D_pos - R_kl; %vector from D-pos to Reflector
        V_DS = S_pos - D_pos; % vector from D-pos to S-pos
        u_in  = V_SR/norm(V_SR); %unit incidence ray
        u_out =  V_RD / norm(V_RD);%unit reflected ray
        n_kl = ((-u_in)+u_out)/sqrt(2+2*((-u_in) * u_out')); %%% normal vector

        e1 = [1, 0, 0];
        e3 = [0, 0, 1]; %3rd column of identity matrix
        beta_kl = asin(n_kl*e3'); 
        gamma_kl =asin((n_kl*e1')/cos(beta_kl));

        V_DR = -V_RD;
        psi_kl_d = atan2(norm(cross(V_DS,V_DR)),dot(V_DS,V_DR)); %% angle between V_DS and V_DR

        G_O_kl = 0; %%%intilization of OIRS GAIN
        if ((psi_kl_d<=psi_FOV_rad) && (psi_kl_d>=0))
            %%% Calculation of G_O_kl
            psi_s_kl = atan2(norm(cross(-V_SR,n_kl)),dot(-V_SR,n_kl));
            c2 = cos(psi_s_kl);

            phi_s_kl = atan2(norm(cross(V_SR,-V_DS)),dot(V_SR,-V_DS));
            c1 = (cos(phi_s_kl))^m_lambert;

            % c3: geometric cosine factor (approx expression preserved from original)
            c3 = (x_kl - x_d)*sin(beta_kl)*cos(gamma_kl)/d_rd + abs(y_kl-y_d)*cos(beta_kl)*cos(gamma_kl)/d_rd + abs(z_kl - z_d)*sin(gamma_kl)/d_rd;

            %c4 is cos(psi_kl_d)
            c4 = cos(psi_kl_d);

            numerator = (rho* eta * (m_lambert+1)*A_PD);
            denominator = 2*pi*((d_sr+d_rd)^2);
            G_O_kl = (numerator/denominator)*c1*c2*c3*c4*T_psi*G_psi;
        else
            G_O_kl = 0;
        end

        h_p = exp(-c_lambda*(d_sr+d_rd)); % attenuation over path

        % Underwater OIRS assisted Channel (deterministic part)
        H_OIRS = H_OIRS + G_O_kl * h_p; 

        count_nonzero = count_nonzero + 1;
    end
end

fprintf('Total deterministic H_OIRS = %.6e\n', H_OIRS);


%% ----- Transmit power sweep & SNR plot (deterministic H_OIRS) -----
Pt_dBm = 0:3:60;
Pt_W   = 1e-3 * 10.^(Pt_dBm/10);          % convert dBm -> W

% Received power using: P_R = Re * Pt * H_OIRS (deterministic)
P_R    = Re .* Pt_W .* H_OIRS;             % vectorized

% SNR and SNR(dB) (deterministic, no turbulence)
SNR    = (P_R).^2 ./ sigma2_i;
SNR_dB = 10*log10(SNR);

figure;
plot(Pt_dBm, SNR_dB, '-o','LineWidth',2,'MarkerFaceColor','auto');
grid on;
xlabel('Transmit Power (dBm)');
ylabel('SNR (dB)');
title('SNR vs Transmit Power using P_R = R_e \cdot P_t \cdot H_{MirrorArray} (no turbulence)');

fprintf('\nPower sweep summary (deterministic H_OIRS):\n');
for i = 1:numel(Pt_dBm)
    fprintf('Pt = %3d dBm -> P_R = %.3e W, SNR = %.2f dB\n', Pt_dBm(i), P_R(i), SNR_dB(i));
end


%% -------- Spectral Efficiency: SE = log2(1 + SNR) (deterministic) --------
SE_log2 = log2(1 + exp(1) / (2* pi) *SNR);   % (user used this earlier)

figure;
plot(Pt_dBm, SE_log2, '-o','LineWidth',2,'MarkerFaceColor','auto');
grid on;
xlabel('Transmit Power (dBm)');
ylabel('Spectral Efficiency (bits/s/Hz)');
title('Spectral Efficiency vs Transmit Power (deterministic)');

fprintf('\nSE summary (log2(1+SNR) approximation):\n');
for i = 1:numel(Pt_dBm)
    fprintf('Pt = %3d dBm -> SE = %.3f bits/s/Hz\n', Pt_dBm(i), SE_log2(i));
end


%% ---------------- BER (no turbulence) -----------
rng(2);                % reproducible
Nbits_no_turb = 1e6;   % bits per Pt (no turbulence)
x_bits = randi([0 1], Nbits_no_turb, 1);   % OOK bitstream (0/1)

H_det = Re * H_OIRS;   % deterministic electrical gain (A/W)

BER_sim = zeros(size(Pt_W));

for ip = 1:numel(Pt_W)
    Pt = Pt_W(ip);           % transmit power in Watts

    % amplitude for symbol '1' at the receiver (current)
    A1 = H_det * Pt;         % (A) - same for every bit since no turbulence

    % transmitted current for each bit (0->0, 1->A1)
    s = A1 .* x_bits;    

    % AWGN noise samples 
    n = sqrt(sigma2_i) .* randn(Nbits_no_turb,1);

    % received
    y = s + n;

    % optimal threshold for OOK in AWGN (equally likely symbols)
    T = 0.5 * A1;

    % detect and compute BER (biterr)
    x_det = y > T;
    [~, BER_sim(ip)] = biterr(x_bits, x_det);
end

% Print a quick table
fprintf('\nPt(dBm)   P_R (W)       SimBER (no turbulence)\n');
for i = 1:numel(Pt_dBm)
    fprintf('%3d       %.3e    %.3e\n', Pt_dBm(i), P_R(i), BER_sim(i));
end

figure;
semilogy(Pt_dBm, BER_sim, 'o-', 'LineWidth', 1.8, 'MarkerFaceColor','auto'); hold on;
grid on;
xlabel('Transmit Power (dBm)');
ylabel('BER');
legend('Simulated BER (no turbulence)');
title('BER vs Transmit Power (OOK, no turbulence)');
ylim([1e-6 1]);


%% ================= MULTI-EGG TURBULENCE + MONTE-CARLO BER =================
% EGG parameter sets from your table (for 3 scintillation indices)

omega_vec  = [0.1665, 0.6238, 0.7210];
lambda_vec = [0.1207, 0.1094, 0.1479];
a_vec      = [0.1559, 0.0111, 0.0121];
b_vec      = [1.5216, 4.4750, 7.4189];
c_vec      = [22.8754, 105.3550, 65.6983];

N_egg = length(omega_vec);   % number of EGG cases

% Monte-Carlo parameters
rng(2);
Nreal = 5000;                 % number of fading realizations
Nbits_per_real = 1e4;

BER_turb_all = zeros(N_egg, numel(Pt_W));
SNR_avg_all_dB = zeros(N_egg, numel(Pt_W));

fprintf('\nRunning MULTI-EGG Monte-Carlo BER simulations...\n');

for ie = 1:N_egg

    omega      = omega_vec(ie);
    lambda_exp = lambda_vec(ie);
    a_gg       = a_vec(ie);
    b_gg       = b_vec(ie);
    c_gg       = c_vec(ie);

    % -------- Generate EGG fading samples --------
    h_alpha = zeros(Nreal,1);
    u = rand(Nreal,1);
    idx_exp = u <= omega;
    idx_gg  = ~idx_exp;

    if any(idx_exp)
        h_alpha(idx_exp) = exprnd(lambda_exp, sum(idx_exp), 1);
    end

    if any(idx_gg)
        Y = gamrnd(a_gg, 1, sum(idx_gg), 1);
        h_alpha(idx_gg) = b_gg .* (Y).^(1./c_gg);
    end

    h_alpha = max(h_alpha, 0);

    % -------- Monte-Carlo BER for this EGG case --------
    x_bits_mat = randi([0 1], Nbits_per_real, Nreal);

    for ip = 1:numel(Pt_W)

        Pt = Pt_W(ip);
        errors_total = 0;
        bits_total = 0;

        for r = 1:Nreal

            h_inst = h_alpha(r);
            H_det_inst = Re * H_OIRS * h_inst;

            A1 = H_det_inst * Pt;
            bits = x_bits_mat(:, r);
            s = A1 .* bits;

            n = sqrt(sigma2_i) .* randn(Nbits_per_real,1);
            y = s + n;

            T = 0.5 * A1;
            x_det = y > T;

            errors_total = errors_total + sum(xor(bits, x_det));
            bits_total = bits_total + Nbits_per_real;
        end

        BER_turb_all(ie, ip) = errors_total / bits_total;
    end

    % -------- Average SNR for this EGG case --------
    SNR_inst = zeros(Nreal, numel(Pt_W));
    for r = 1:Nreal
        for ip = 1:numel(Pt_W)
            A1 = Re * Pt_W(ip) * H_OIRS * h_alpha(r);
            SNR_inst(r, ip) = (A1.^2) / sigma2_i;
        end
    end

    SNR_avg_all_dB(ie,:) = 10*log10(mean(SNR_inst,1));

    fprintf('EGG Case %d completed.\n', ie);

end

%% ================= FINAL BER vs SNR COMPARISON PLOT =================
figure;

% --- No turbulence ---
semilogy(SNR_dB, BER_sim, 'ko-', 'LineWidth', 2, 'MarkerFaceColor','auto'); 
hold on;

% --- EGG turbulence cases ---
markers = {'s','^','d'};
for ie = 1:N_egg
    semilogy(SNR_avg_all_dB(ie,:), BER_turb_all(ie,:), ...
        [markers{ie} '-'], 'LineWidth', 2, 'MarkerFaceColor','auto');
end

grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('No Turbulence', ...
       'EGG – Weak Turbulence', ...
       'EGG – Moderate Turbulence', ...
       'EGG – Strong Turbulence', ...
       'Location','southwest');
title('Mirror Array: BER vs SNR for without Turbulence and with EGG Turbulence');
ylim([1e-6 1]);

fprintf('\nMULTI-EGG Simulation Completed.\n');

