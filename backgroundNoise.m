
function bg = backgroundNoise(t, lat, lon)
% Estimate Sun elevation angle without Aerospace Toolbox
% Inputs: t (datetime, UTC), lat/lon in deg
% Output: background count rate (cps)

% --- convert to Julian date ---
jd = juliandate(t);
d = jd - 2451545.0;  % days since J2000

% --- mean longitude and anomaly ---
g = deg2rad(mod(357.529 + 0.98560028*d, 360));
q = deg2rad(mod(280.459 + 0.98564736*d, 360));
L = q + deg2rad(1.915)*sin(g) + deg2rad(0.020)*sin(2*g);

% --- obliquity ---
eps = deg2rad(23.439 - 0.00000036*d);

% --- RA and declination ---
RA = atan2(cos(eps)*sin(L), cos(L));
dec = asin(sin(eps)*sin(L));

% --- local sidereal time ---
T = (jd - floor(jd)) * 24; % UT hours
theta = deg2rad(mod(280.16 + 360.9856235*(d + T/24), 360));
H = theta + deg2rad(lon) - RA; % hour angle

% --- Sun elevation ---
latr = deg2rad(lat);
el = asin(sin(latr)*sin(dec) + cos(latr)*cos(dec)*cos(H));

% --- background model ---
if el > 0
    bg = 800; % cps daytime
else
    % moonlight proxy
    phase = mod(day(t),29)/29;
    if phase > 0.5
        bg = 100; % cps bright night (moon)
    else
        bg = 10;  % cps dark night
    end
end
end
