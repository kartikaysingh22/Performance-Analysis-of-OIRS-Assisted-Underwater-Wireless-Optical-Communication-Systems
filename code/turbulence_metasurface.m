% ============================================================
%   OIRS-assisted UWOC Simulation
%   Underwater Wireless Communication System"
%    METASURFACES TURBULENCE EGG
% ============================================================

clear; clc; close all;

%% -------------------- 1. System Parameters --------------------
% (Use exact paper Table 2 values here if available)

A_o   = 1;          % Total OIRS area (m^2)
n_elem2 = 169;      % Total number of mirror elements (n_m^2)
l_m   = 0.08;       % Mirror element length (m)
w_m   = 0.08;       % Mirror element width (m)
y_s   = 1;         % Tx depth from OIRS (m)
y_d   = 1;         % Rx depth from OIRS (m)
n_p = 13;           %%% no of patches 13 x 13

% --- EGG turbulence parameters (1st values) ---
omega = 0.1665;        % mixture weight
lambda_exp = 0.1207;   % parameter
a_gg = 0.1559;         % gamma param a
b_gg = 1.5216;         % gamma param b
c_gg = 22.8754;         %  gamma param c

% --- Physical / receiver constants ---
A_PD  = 0.0001;       % photodiode area (m^2) (same as Sipm for now)
psi_FOV_deg = 90;   % receiver FoV (deg)
psi_FOV_rad = deg2rad(psi_FOV_deg); %conversion into radians
mu = 1.33;          % refractive index (water)
rho_MS = 0.8;          % mirror reflection coefficient
eta = 0.95;            % LED efficiency
T_psi = 1;          % optical filter gain
G_psi = 1;
c_lambda = 0.056;      % attenuation coefficient (m^-1)
Re = 0.5;           % Re of photodiode
sigma2_i = 1e-12;   % total noise variance (A^2)
c_water = 2.25e8;   % light speed in water (m/s)

% --- LED Lambertian source parameters ---
theta_half_deg = 60;  %% taken from Vivek Sir's and Anand Singh's paper
m_lambert = 2; % m_lambert= -log(2)/log(cosd(theta_half_deg)) fomula


%% -------------------- 2. Geometry --------------------

d_sd = 4; % distance between source and destination
% Tx-Rx separation (m)
x_s = 0; %%% x-cordinate of source led
x_d = x_s; % x cordinate of Rx
z_d = d_sd/2; % z cordinate of Rx
S_pos = [x_s, y_s, -d_sd/2];    % Tx position
D_pos = [x_s, y_d, z_d];    % Rx position



%% -------------------- 4. Compute OIRS Channel Gain --------------------
fprintf('Computing per-element OIRS gains...\n');
H_OIRS = 0;
count_nonzero = 0;
total_length = n_p * l_m;
total_width  = n_p * w_m;
x_0 = -0.56; %% k will update x
z_0 = -0.56; %% l will update z
for k = 1:n_p
    x_kl = x_0 + 0.08 * k;
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
        n_kl = [0,1,0] % normal of a patch remains same for each patch
        
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
        
        else
            G_O_kl = 0; %%% other cases
        end

        h_p = exp(-c_lambda*(d_sr+d_rd)); %%% Attentuation Constant
        

        %%% Underwater OIRS assited Channel
        H_OIRS = H_OIRS + G_O_kl * h_p;
        
        count_nonzero = count_nonzero + 1;
    end
    z_kl = -0.56;
end

fprintf('Total H_OIRS = %.6e\n', H_OIRS);

%% ----- Transmit power sweep & SNR plot -----
% Power vector (dBm) and Watts
Pt_dBm = 0:3:60;
Pt_W   = 1e-3 * 10.^(Pt_dBm/10);          % convert dBm -> W

% Received power using: P_R = Re * Pt * H_OIRS
P_R    = Re .* Pt_W .* H_OIRS;             % vectorized

% SNR and SNR(dB)
SNR    = (P_R).^2 ./ sigma2_i;
SNR_dB = 10*log10(SNR);

% Plot
figure;
plot(Pt_dBm, SNR_dB, '-o','LineWidth',2,'MarkerFaceColor','auto');
grid on;
xlabel('Transmit Power (dBm)');
ylabel('SNR (dB)');
title('SNR vs Transmit Power using P_R = R_e \cdot P_t \cdot H_{Metasurface}');

% (Optional) print a quick summary
fprintf('\nPower sweep summary:\n');
for i = 1:numel(Pt_dBm)
    fprintf('Pt = %3d dBm -> P_R = %.3e W, SNR = %.2f dB\n', ...
        Pt_dBm(i), P_R(i), SNR_dB(i));
end



%% -------------------- 5. Spectral Efficiency (Shannon formula) --------------------
% SE = log2(1 + SNR)

SE_log2 = log2(1 + exp(1) / (2* pi) * SNR);   % Spectral efficiency in bits/s/Hz

% Plot SE vs transmit power
figure;
plot(Pt_dBm, SE_log2, '-o','LineWidth',2,'MarkerFaceColor','auto');
grid on;
xlabel('Transmit Power (dBm)');
ylabel('Spectral Efficiency (bits/s/Hz)');
title('Spectral Efficiency vs Transmit Power');

% Display quick summary in Command Window
fprintf('\nSpectral Efficiency summary (log2(1 + SNR)):\n');
for i = 1:numel(Pt_dBm)
    fprintf('Pt = %3d dBm -> SE = %.4f bits/s/Hz\n', Pt_dBm(i), SE_log2(i));
end






%% ---------------- BER (no turbulence) - Simulated  -----------


rng(2);                % reproducible
Nbits = 1e6;           % bits per Pt 
x_bits = randi([0 1], Nbits, 1);   % OOK bitstream (0/1)

H_det = Re * H_OIRS;   % deterministic electrical gain (A/W) already includes attenuation

BER_sim = zeros(size(Pt_W));

for ip = 1:numel(Pt_W)
    Pt = Pt_W(ip);           % transmit power in Watts

    % amplitude for symbol '1' at the receiver (current)
    A1 = H_det * Pt;         % (A) - same for every bit since no turbulence

    % transmitted current for each bit (0->0, 1->A1)
    s = A1 .* x_bits;        % Nbits x 1

    % AWGN noise samples 
    n = sqrt(sigma2_i) .* randn(Nbits,1);

    % received
    y = s + n;

    % optimal threshold for OOK in AWGN (equally likely symbols)
    T = 0.5 * A1;

    % detect and compute BER (biterr)
    x_det = y > T;
    [~, BER_sim(ip)] = biterr(x_bits, x_det);

end

% Print a quick table
fprintf('\nPt(dBm)   P_R (W)       SimBER\n');
for i = 1:numel(Pt_dBm)
    fprintf('%3d       %.3e    %.3e    %.3e\n', Pt_dBm(i), P_R(i), BER_sim(i), 'n');
end

% Plot simulated and theoretical BER vs Pt (log scale)
figure;
semilogy(Pt_dBm, BER_sim, 'o-', 'LineWidth', 1.8, 'MarkerFaceColor','auto'); hold on;
grid on;
xlabel('Transmit Power (dBm)');
ylabel('BER');
legend('Simulated BER (no turbulence)');
title('BER vs Transmit Power (OOK, no turbulence)');
ylim([1e-6 1]); 


%% ---------------- EGG turbulence and Monte Carlo BER ----------------
% Monte-Carlo parameters (changeable)
rng(2);                     % reproducible
Nreal = 10000;               % number of channel realizations/simulations (Monte-Carlo)
Nbits_per_real = 1e4;       % size of sample bitstream

% EGG samples
h_alpha = zeros(Nreal,1);

% choosing when to use exp and when to use gamma distribution
u = rand(Nreal,1);
idx_exp = u <= omega;
idx_gg  = ~idx_exp;

% Exponential samples (scale = lambda_exp)
if any(idx_exp)
    h_alpha(idx_exp) = exprnd(lambda_exp, sum(idx_exp), 1);
end

% Generalized-Gamma samples: If Y ~ Gamma(a,1) then X = b * Y^(1/c)
if any(idx_gg)
    Y = gamrnd(a_gg, 1, sum(idx_gg), 1);   
    h_alpha(idx_gg) = b_gg .* (Y).^(1./c_gg);
end

% Floor at zero to prevent from h_alpha going negative
h_alpha = max(h_alpha, 0); 

% Calculating simulated scintillation index (SI)
mean_h = mean(h_alpha);
var_h  = var(h_alpha);
SI = var_h / (mean_h^2); % Normalized Variance
fprintf('EGG fading: mean(h)=%.4e, var(h)=%.4e, SI=%.4e (Nreal=%d)\n', mean_h, var_h, SI, Nreal);

% Showing the PDFof the simulated EGG
figure;
histogram(h_alpha,100,'Normalization','pdf'); hold on;
xlabel('h (fading gain)');
ylabel('PDF');
title('Sampled EGG PDF (histogram)');
% overlay KDE(kernel Density Estimate)
try
    [f,xi] = ksdensity(h_alpha);
    plot(xi,f,'LineWidth',1.5);
    legend('Histogram','KDE');
catch
    legend('Histogram');
end

%% Monte-Carlo BER under EGG turbulence
BER_turb = zeros(size(Pt_W));
rng(2);
% Pre-generating bitstreams for each realization
% matrix
% each column is a sample bitstream
x_bits_mat = randi([0 1], Nbits_per_real, Nreal);

fprintf('Running BER Monte-Carlo under EGG turbulence...\n');
for ip = 1:numel(Pt_W)
    Pt = Pt_W(ip);
    errors_total = 0;
    bits_total = 0;
    %%%% running Nreal sims to get an average value of BER 
    %%%% (Monte Carlo Sim technique)
    for r = 1:Nreal
        h_inst = h_alpha(r);                  % fading sample for this realization
        H_det_inst = Re * H_OIRS * h_inst;  % instantaneous electrical gain (A/W)

        % amplitude for symbol '1' at the receiver (current)
        A1 = H_det_inst * Pt;    % (A)

        % bits for this realization
        bits = x_bits_mat(:, r);
        s = A1 .* bits;          % transmitted current for each bit (0->0, 1->A1)

        % AWGN noise samples
        n = sqrt(sigma2_i) .* randn(Nbits_per_real,1);

        % Received
        y = s + n;

        % Threshold (optimal for equally likely OOK)
        T = 0.5 * A1;

        % Detection and errors
        x_det = y > T;
        errors_total = errors_total + sum(xor(bits, x_det));
        bits_total = bits_total + Nbits_per_real;
    end
    BER_turb(ip) = errors_total / bits_total;
    fprintf('Pt = %2d dBm -> BER_turb = %.3e\n', Pt_dBm(ip), BER_turb(ip));
end

% Plot BER (turbulent) and compare to no-turbulence curve
figure;
semilogy(Pt_dBm, BER_turb, 's-','LineWidth',1.8,'MarkerFaceColor','auto'); hold on;
semilogy(Pt_dBm, BER_sim, 'o-','LineWidth',1.2,'MarkerFaceColor','auto');
legend('BER (EGG turbulence)','BER (no turbulence)');
xlabel('Transmit Power (dBm)');
ylabel('BER');
grid on;
title('BER vs Transmit Power (OOK) with EGG turbulence');
ylim([1e-6 1]);

%% -------------- Average SNR (with EGG turbulence) --------------
rng(2);
SNR_inst = zeros(Nreal, numel(Pt_W));
for r = 1:Nreal
    for ip = 1:numel(Pt_W)
        A1 = Re * Pt_W(ip) * H_OIRS * h_alpha(r);
        SNR_inst(r, ip) = (A1.^2) / sigma2_i;
    end
end
SNR_avg_dB = 10*log10(mean(SNR_inst, 1));

figure;
plot(Pt_dBm, SNR_avg_dB, '-o','LineWidth',1.6);
xlabel('Transmit Power (dBm)');
ylabel('Average SNR (dB)');
grid on;
title('Average SNR across EGG realizations');




%% ------------------BER vs SNR (No Turbulence & EGG Turbulence)------------------------- %%

figure;
semilogy(SNR_avg_dB, BER_turb, 'rs-', 'LineWidth', 1.8, 'MarkerFaceColor','auto');
hold on;
semilogy(SNR_dB, BER_sim, 'bo-', 'LineWidth', 1.8, 'MarkerFaceColor','auto'); 
hold on;

grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('BER (EGG Turbulence)', 'BER (No Turbulence)','Location','southwest');
title('BER vs SNR (OOK) for No Turbulence and EGG Turbulence');

ylim([1e-6 1]);
xlim([min([SNR_dB, SNR_avg_dB]) max([SNR_dB, SNR_avg_dB])]);

fprintf('\nSim Completed.\n');


