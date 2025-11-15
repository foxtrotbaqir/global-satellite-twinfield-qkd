%% TF-QKD MULTIHOP ROUTING WITH MICIUS-STYLE CONSTELLATION
% With: asymmetric decoys, realistic losses, finite-key effects,
% and Sun/Moon-driven background noise
%
% Live CSV columns: time, path1, skr1_bps, path2, skr2_bps, path3, skr3_bps, skr_xor_bps

clear; clc; close all;

%% PARAMETERS
n_orbits       = 25;                % Orbital planes
sats_per_orbit = 44;                % Satellites per plane
altitude_km    = 500;               % Altitude [km]
Re             = 6371e3;            % Earth radius [m]
h              = altitude_km * 1e3; % Orbit altitude [m]
a              = Re + h;            % Semi-major axis [m]
ecc            = 0.0011;            % Slight eccentricity
inc            = 97.369;            % Inclination from Micius TLE
raan0          = 229.5573;          % Reference RAAN
raan_gap       = 360 / n_orbits;    % RAAN spacing [deg]
argPerigee     = 112.8982;          % Argument of perigee
total_sats     = n_orbits * sats_per_orbit;

% Timing
startTime  = datetime(2025,7,31,12,0,0,'TimeZone','UTC');
stopTime   = startTime + days(30);
sampleTime = 1;                      % propagation resolution
stepSec    = 1;                     % dataset cadence
timeSteps  = startTime:seconds(stepSec):stopTime;

sc_path = satelliteScenario(startTime, stopTime, sampleTime);

% Ground stations
gs1 = groundStation(sc_path, 48.1351, 11.5820, 'Name', 'Munich');
gs2 = groundStation(sc_path, 40.7128, -74.0060, 'Name', 'NewYork');
gs1.MinElevationAngle = 20; gs2.MinElevationAngle = 20; %tighter elevation mask for better efficiency

%% BUILD CONSTELLATION
satNames = strings(1,total_sats);
satObjs(total_sats) = satellite(sc_path,a,ecc,inc,raan0,argPerigee,0,'Name','Temp');
for i = 0:n_orbits-1
    raan = mod(raan0 + i*raan_gap,360);
    for j = 0:sats_per_orbit-1
        idx = i*sats_per_orbit + j + 1;
        ta  = j*(360/sats_per_orbit);
        nm  = sprintf('Sat-%02d-%02d',i+1,j+1);
        satObjs(idx)  = satellite(sc_path,a,ecc,inc,raan,argPerigee,ta,'Name',nm);
        satNames(idx) = nm;
    end
end

%% CSV
CSV_FILE  = "tfqkd_paths_xor_dataset.csv";
if ~isfile(CSV_FILE)
    fid=fopen(CSV_FILE,'w');
    fprintf(fid,'time,path1,skr1_bps,path2,skr2_bps,path3,skr3_bps,skr_xor_bps\n');
    fclose(fid);
end

%% MAIN LOOP
SKRaccum = zeros(1,numel(timeSteps));
for ti = 1:numel(timeSteps)
    t = timeSteps(ti);
    fprintf("\nTime %s:\n",string(t));

    % Node graph
    nodeNames = ['Munich', satNames, 'NewYork'];
    G = graph; G = addnode(G,nodeNames);

    gs1Pos = myLla2ecef([gs1.Latitude, gs1.Longitude, 0]);
    gs2Pos = myLla2ecef([gs2.Latitude, gs2.Longitude, 0]);

    satPos = zeros(3,total_sats);
    for k=1:total_sats
        satPos(:,k) = states(satObjs(k),t,'CoordinateFrame','ecef');
    end

    % GS–Sat
    for k=1:total_sats
        [d1,ang1] = rangeangle(satPos(:,k),gs1Pos); el1=ang1(2);
        [d2,ang2] = rangeangle(satPos(:,k),gs2Pos); el2=ang2(2);
        if d1<1200e3 && el1>20, G=addedge(G,'Munich',satNames(k),d1/1e3); end
        if d2<1200e3 && el2>20, G=addedge(G,satNames(k),'NewYork',d2/1e3); end
    end
    % ISLs
    for iSat=1:total_sats-1
        for jSat=iSat+1:total_sats
            d=norm(satPos(:,iSat)-satPos(:,jSat));
            if d<1500e3, G=addedge(G,satNames(iSat),satNames(jSat),d/1e3); end
        end
    end

    % Node-disjoint shortest paths
    paths={}; Gc=G;
    for p=1:3
        [sp,~]=shortestpath(Gc,'Munich','NewYork');
        if isempty(sp), break; end
        paths{end+1}=sp; %#ok
        Gc=rmnode(Gc,sp(2:end-1));
    end
    fprintf("  %d disjoint paths found.\n",numel(paths));

    % Per-path SKRs
    pathStrs=repmat("(none)",1,3); perPathSKR=zeros(1,3);
    for pi=1:min(3,numel(paths))
        namesPi=paths{pi};
        perPathSKR(pi)=computePathSKR_asym(namesPi,t,satObjs,gs1,gs2,satNames);
        pathStrs(pi)=strjoin(string(namesPi),'->');
        fprintf("  Path %d SKR = %.2f bps\n",pi,perPathSKR(pi));
    end

    % XOR
    nz=perPathSKR(perPathSKR>0);
    skr_xor=0; if ~isempty(nz), skr_xor=min(nz); end
    fprintf("  XOR SKR = %.2f bps\n",skr_xor);
    SKRaccum(ti)=skr_xor;

    % Append CSV
    rowT=table(string(t),pathStrs(1),perPathSKR(1),pathStrs(2),perPathSKR(2), ...
               pathStrs(3),perPathSKR(3),skr_xor, ...
               'VariableNames',{'time','path1','skr1_bps','path2','skr2_bps','path3','skr3_bps','skr_xor_bps'});
    writetable(rowT,CSV_FILE,'WriteMode','Append','WriteVariableNames',false);
end

% Plot
figure; bar(timeSteps,SKRaccum,'FaceColor',[0.2 0.5 0.8],'EdgeColor','none');
xlabel("Time (UTC)"); ylabel("XOR SKR [bps]");
title("TF-QKD SKR Over Time with Finite-Key + Background Noise"); grid on;
saveas(gcf, 'tfqkd_skr_barplot.png');


