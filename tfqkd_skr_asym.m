function R=tfqkd_skr_asym(etaA,etaB,t,gs)
% Constants
Clock_rate=1e9; d=0.5; M=16; eta_sync=0.95;
f=1.15; e_opt=0.03; E_M=0.01275;
P_dc=10/Clock_rate; % intrinsic dark count

% Background noise model (Sun/Moon)
bg=backgroundNoise(t,gs.Latitude,gs.Longitude);
P_bg=bg/Clock_rate;
P_tot=P_dc+P_bg;

% Decoy intensities (optimized scaling)
muB=0.5; nuB=0.1; wB=1e-4;
muA=muB*(etaB/etaA); nuA=nuB*(etaB/etaA); wA=wB;

% Gains
Qmu=1-(1-P_tot)^2*exp(-(muA*etaA+muB*etaB));
Qnu=1-(1-P_tot)^2*exp(-(nuA*etaA+nuB*etaB));
Qw =1-(1-P_tot)^2*exp(-(wA*etaA +wB*etaB));

% Errors
Emu=0.5+(1/(2*Qmu))*(1-P_tot)*(exp(-(muA*etaA+muB*etaB)*(1-(e_opt+E_M)))- ...
                                     exp(-(muA*etaA+muB*etaB)*(e_opt+E_M)));

% Yields
y0=(nuB*Qw*exp(wB)-wB*Qnu*exp(nuB))/(nuB-wB);
y1=((muB^2)*Qnu*exp(nuB)-(muB^2)*Qw*exp(wB)-(nuB^2-wB^2)*(Qmu*exp(muB)-y0)) ...
    /(muB*(muB*nuB-muB*wB-nuB^2+wB^2));
e1=(Emu*Qmu*exp(muB)-0.5*y0)/(y1*muB);
Q1=exp(-muB)*muB*y1;

% Binary entropy
h=@(x) -x.*log2(x)-(1-x).*log2(1-x);

% Finite-key correction (block size N)
N=1e9; % pulses/block
delta=5*sqrt(log(2/1e-10)/N); % security parameter
R_QKD=Q1*(1-h(e1-delta))-f*Qmu*h(Emu+delta);

R=max(0,eta_sync*(d*Clock_rate/M)*R_QKD);
end
