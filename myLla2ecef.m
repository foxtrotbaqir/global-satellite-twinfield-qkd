function ecef = myLla2ecef(lla)
% Convert [lat, lon, alt] (deg, deg, m) to ECEF [x,y,z] (m)
% No toolboxes required (WGS84)

lat = deg2rad(lla(1));
lon = deg2rad(lla(2));
alt = lla(3);

% WGS84 constants
a = 6378137.0;          % semi-major axis [m]
f = 1/298.257223563;    % flattening
e2 = f*(2-f);           % eccentricity squared

% Prime vertical radius of curvature
N = a / sqrt(1 - e2*sin(lat)^2);

% ECEF coordinates
x = (N + alt) * cos(lat) * cos(lon);
y = (N + alt) * cos(lat) * sin(lon);
z = ((1 - e2)*N + alt) * sin(lat);

ecef = [x; y; z];
end