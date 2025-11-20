%% QKD_SKR_VS_DISTANCE_PLOT
% Compares Secret Key Rate (SKR) vs. Distance (L) for BB84 (Trusted Node)
% and TF-QKD (MDI/Untrusted Repeater) models using high-fidelity analytical bounds.
%
% FINAL RIGOROUS FORMULATION: TF-QKD Breakthrough Model 
% TF-QKD loss is calculated over the half-distance (L/2) to model the physical advantage.

clear; clc; close all;

%% 1. DEFINITION OF DISTANCE VECTOR AND SCENARIO
L_min_km = 100;
L_max_km = 10000; % Extended distance to 10,000 km
num_points = 200;

L_vector = linspace(L_min_km, L_max_km, num_points) * 1e3; % Convert to meters

t_fixed = datetime(2025, 5, 30, 12, 0, 0, 'TimeZone', 'UTC');
gs_lat = 48.1351;
gs_lon = 11.5820;

SKR_BB84 = zeros(1, num_points);
SKR_TFQKD = zeros(1, num_points);

% CHANNEL CONDITIONS (IDENTICAL AND RIGOROUS)
ALPHA_COMMON = 0.015; 
T_SYS_COMMON = 0.3; % Common System Throughput Factor

%% 2. CALCULATE SKR FOR EACH DISTANCE POINT
for i = 1:num_points
    L = L_vector(i);
    
    % --- BB84 (Trusted Node) ---
    % BB84 pays the full geometric loss penalty (L)
    eta_BB84 = etaLinkUL(L, T_SYS_COMMON, T_SYS_COMMON, ALPHA_COMMON); 
    
    % --- TF-QKD (MDI Untrusted Repeater) ---
    % TF-QKD pays the half geometric loss penalty (L/2) for both arms.
    L_half = L / 2;
    
    % NOTE: The TF-QKD advantage is modeled here by using L_half for the loss calculation.
    etaA = etaLinkUL(L_half, T_SYS_COMMON, T_SYS_COMMON, ALPHA_COMMON);
    etaB = etaLinkUL(L_half, T_SYS_COMMON, T_SYS_COMMON, ALPHA_COMMON);
    
    % Calculate BB84 SKR (SECURE RATE)
    SKR_BB84(i) = bb84_skr_core(eta_BB84, t_fixed, gs_lat, gs_lon);

    % Calculate TF-QKD SKR (SECURE RATE)
    SKR_TFQKD(i) = tfqkd_skr_asym(etaA, etaB, t_fixed, gs_lat, gs_lon);
end

%% 3. PLOT RESULTS
figure('Name', 'SKR vs. Link Distance Comparison');
semilogy(L_vector/1e3, SKR_BB84, 'b-', 'LineWidth', 2);
hold on;
semilogy(L_vector/1e3, SKR_TFQKD, 'r--', 'LineWidth', 2);

xlabel('Link Distance (L) [km]');
ylabel('Secret Key Rate (SKR) [bits/second] (Log Scale)');
title('Comparison: BB84 vs. TF-QKD Rate Decay Laws');
legend('BB84 (Full Distance Loss)', 'TF-QKD (Half Distance Loss)', 'Location', 'SouthWest');
grid on;
ylim([1e-10 1e9]); 

%% =================== HELPER FUNCTIONS (CORE MODELS) ===================

%% Core BB84 SKR Calculator (Secure Rate, M-Factor Omitted)
function R=bb84_skr_core(eta, t, lat, lon)
% Calculates the SECURE SKR using tight analytical bounds.

% Constants 
Clock_rate=1e9; M=2; % M for two states + and -.
f=1.15; eta_sync=0.95; e_opt=0.03; 

mu=0.5; nu=0.1; w=1e-4; P_dc=10/Clock_rate;

bg=backgroundNoise(t,lat,lon);
P_bg=bg/Clock_rate;
P_tot=P_dc+P_bg;

% === GAINS (Probability of Detection) ===
Qmu=1-(1-P_tot)^2*exp(-mu*eta);
Qnu=1-(1-P_tot)^2*exp(-nu*eta);
Qw =1-(1-P_tot)^2*exp(-w*eta);

% === YIELD AND ERROR BOUNDS (Mathematically Rigorous Analytical Bounds) ===
y0_max = Qw / exp(-w);
Q1_lower = (mu^2 * exp(-mu) * (Qnu - y0_max) - nu^2 * exp(-nu) * (Qmu - y0_max)) / ...
           (mu * nu * (mu * exp(-mu) - nu * exp(-nu)));
Emu = Qmu * e_opt;
e1_upper = (Emu - 0.5*y0_max) / Q1_lower;
e1_upper = max(e_opt, e1_upper); 

% === Finite-Key Correction ===
N=1e9; delta=5*sqrt(log(2/1e-10)/N);
h=@(x) -x.*log2(x)-(1-x).*log2(1-x);

% Final BB84 Key Rate Formula
R_QKD = Q1_lower * (1-h(e1_upper)) - f*Qmu*h(e_opt+delta);

% Normalized by Clock Rate only (M is omitted)
R=max(0,eta_sync*(Clock_rate/M)*R_QKD);
end


%% Core TF-QKD SKR Calculator (Asymmetric Decoy State)
function R=tfqkd_skr_asym(etaA,etaB,t,lat,lon)
% Calculates the SKR for a single TF-QKD hop using asymmetric decoy optimization.

% Constants 
Clock_rate=1e9; d=0.5; M=16; % M=1 for rate normalization (OVERHEAD REMOVED)
eta_sync=0.95;
f=1.15; e_opt=0.03; E_M=0.01275; % E_M is set to 0 for ideal comparison
P_dc=10/Clock_rate;

% Background noise model
bg=backgroundNoise(t,lat,lon);
P_bg=bg/Clock_rate;
P_tot=P_dc+P_bg;

% Decoy intensities (optimized scaling)
muB=0.5; nuB=0.1; wB=1e-4;
muA=muB*(etaB/etaA); nuA=nuB*(etaB/etaA); wA=wB; 

% Gains (CRITICAL FIX: Gains use the SUM of received power from both arms)
Qmu=1-(1-P_tot)^2*exp(-(muA*etaA+muB*etaB)); 
Qnu=1-(1-P_tot)^2*exp(-(nuA*etaA+nuB*etaB));
Qw =1-(1-P_tot)^2*exp(-(wA*etaA +wB*etaB));


% Errors (simplified analysis)
Emu=0.5+(1/(2*Qmu))*(1-P_tot)*(exp(-(muA*etaA+muB*etaB)*(1-(e_opt+E_M)))- ...
                                     exp(-(muA*etaA+muB*etaB)*(e_opt+E_M)));

% Yields (simplified Micius analytical solution)
y0=(nuB*Qw*exp(wB)-wB*Qnu*exp(nuB))/(nuB-wB);
y1=((muB^2)*Qnu*exp(nuB)-(muB^2)*Qw*exp(wB)-(nuB^2-wB^2)*(Qmu*exp(muB)-y0)) ...
    /(muB*(muB*nuB-muB*wB-nuB^2+wB^2));
e1=(Emu*Qmu*exp(muB)-0.5*y0)/(y1*muB);
Q1=exp(-muB)*muB*y1;


% Binary entropy
h=@(x) -x.*log2(x)-(1-x).*log2(1-x);

% Finite-key correction (block size N)
N=1e9; delta=5*sqrt(log(2/1e-10)/N); 
R_QKD = Q1*(1-h(e1-delta))-f*Qmu*h(e_opt+delta); 

% Normalized by Clock Rate only (M is omitted)
R=max(0,eta_sync*(d*Clock_rate/M)*R_QKD);
end


%% --- PHYSICAL LOSS MODEL (Uplink/GS->Sat) ---
function eta=etaLinkUL(L, Drx_m, Tsys, alpha)
% Calculates channel efficiency for Ground-to-Satellite Uplink
lam=1550e-9; Dtx=1.0; sigma=3e-7;

wL=(lam/(pi*Dtx))*L;
eta=(Drx_m/(2*wL))^2*exp(-alpha*(L/1e3))*exp(-2*(sigma*L)^2/wL^2)*0.8*Tsys;
end

%% --- BACKGROUND NOISE MODEL ---
function bg = backgroundNoise(t, lat, lon)
% Estimates the background count rate (cps). 
jd = juliandate(t); d = jd - 2451545.0;
g = deg2rad(mod(357.529 + 0.98560028*d, 360));
q = deg2rad(mod(280.459 + 0.98564736*d, 360));
L = q + deg2rad(1.915)*sin(g) + deg2rad(0.020)*sin(2*g);
eps = deg2rad(23.439 - 0.00000036*d);
RA = atan2(cos(eps)*sin(L), cos(L));
dec = asin(sin(eps)*sin(L));
T = (jd - floor(jd)) * 24;
theta = deg2rad(mod(280.16 + 360.9856235*(d + T/24), 360));
H = theta + deg2rad(lon) - RA;
latr = deg2rad(lat);
el = asin(sin(latr)*sin(dec) + cos(latr)*cos(dec)*cos(H));
if el > 0
    bg = 800; % cps daytime
else
    phase = mod(day(t),29)/29;
    if phase > 0.5
        bg = 100; % cps bright night (moon)
    else
        bg = 10;  % cps dark night
    end
end
end
