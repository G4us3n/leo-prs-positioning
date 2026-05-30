function ctx = stage01_constellation(cfg, ctx)
%STAGE01_CONSTELLATION  Build the satellite scenario, constellation and ground stations.
%
%   Stage of the LEO-PRS positioning pipeline. Reads the shared context
%   struct CTX (and the CFG parameters), runs one logical section of the
%   simulation, and returns CTX with new fields added. See config.m for
%   parameters and main.m for the run order.

% Pull configuration into local names matching the original script
zenithTime    = cfg.zenithTime;
a             = cfg.altitude_m;
ecc           = cfg.eccentricity;
argp          = cfg.argPerigee;
centerLat     = cfg.centerLat;
centerLon     = cfg.centerLon;
nPlanes       = cfg.nPlanes;
nSatsPerPlane = cfg.nSatsPerPlane;
base_inc      = cfg.inclination_deg;
plane_raan_offsets = cfg.plane_raan_offsets;
train_nu_offsets   = cfg.train_nu_offsets;
cityNetwork   = cfg.cityNetwork;
nSats         = nPlanes * nSatsPerPlane;
scenarioStart = zenithTime - minutes(cfg.windowMinutes);
scenarioStop  = zenithTime + minutes(cfg.windowMinutes);
sc = satelliteScenario(scenarioStart, scenarioStop, cfg.sampleTime);
numCities = size(cityNetwork, 1);

% ----- original section body (unchanged physics) ---------------------
% ---- Generate satellite names ----
% Format: P<plane>-S<sat>, e.g.  P1-S1, P2-S3, P4-S2
satNames = cell(1, nSats);
for p = 1:nPlanes
    for s = 1:nSatsPerPlane
        idx = (p-1)*nSatsPerPlane + s;
        satNames{idx} = sprintf('P%d-S%d', p, s);
    end
end

fprintf('--> CONSTELLATION: %d planes × %d sats = %d total\n', nPlanes, nSatsPerPlane, nSats);
fprintf('    Plane RAAN offsets: [%s] deg\n', num2str(plane_raan_offsets, '%+.1f '));
fprintf('    Train nu offsets:   [%s] deg\n', num2str(train_nu_offsets, '%+.1f '));

% ---- Orbital-mechanics constants ----
mu_earth   = 3.986004418e14;           % Earth's gravitational parameter [m³/s²]
n_mean     = sqrt(mu_earth / a^3);     % Mean motion [rad/s]:
                                        %   angular velocity of the orbit.
                                        %   For a = 6 871 km → n ≈ 0.00112 rad/s
                                        %   → orbital period ≈ 94 min.
timeOffset = seconds(zenithTime - scenarioStart);  % 180 s (3 min)
delta_nu   = rad2deg(n_mean * timeOffset);
    % delta_nu: how many degrees the satellite travels in 3 minutes.
    % ≈ 11.5° for our orbit.

% WGS84: the standard mathematical model of Earth's shape (ellipsoid).
wgs84 = wgs84Ellipsoid('meter');

% ---- Compute the base RAAN so that satellites pass over Milan ----
%
% PROBLEM:  We want all satellites overhead Milan at zenithTime.
%           We must work backwards to find the correct starting Keplerian
%           elements at scenarioStart (3 minutes earlier).
%
% STEP 1:  Find the true anomaly (nu) at which the orbit's latitude
%          matches Milan's latitude.
%          From spherical trig:  sin(lat) = sin(inc) × sin(nu)
%          →  nu_target = asin(sin(lat) / sin(inc))
inc       = base_inc;
nu_target = asind(sind(centerLat) / sind(inc));   % ≈ 56.7° for Milan

% STEP 2:  Rewind by delta_nu to get the starting true anomaly.
nu_start_exact = mod(nu_target - delta_nu, 360);

% STEP 3:  Create a temporary satellite with RAAN = 0° and propagate it
%          to zenithTime.  Compare its sub-satellite longitude with Milan
%          to compute the required RAAN correction.
sat_test = satellite(sc, a, ecc, inc, 0, argp, nu_start_exact);
[pos_test, ~] = states(sat_test, zenithTime, "CoordinateFrame", "ecef");
[~, testLon, ~] = ecef2geodetic(wgs84, pos_test(1), pos_test(2), pos_test(3));
raan_base = mod(centerLon - testLon, 360);   % RAAN shift to hit Milan's longitude
delete(sat_test);                             % Clean up the test satellite

% ---- Create all 12 satellites ----
% MATLAB quirk: satellite() returns a Satellite object.  You cannot index
% into an uninitialised double array with it (e.g. sats = []; sats(2) = satellite(...))
% would fail.  So we create the first satellite separately to establish the
% correct object-array type, then append the rest inside the loop.
final_raan_1 = raan_base + plane_raan_offsets(1);
final_nu_1   = nu_start_exact + train_nu_offsets(1);
sats = satellite(sc, a, ecc, base_inc, final_raan_1, argp, final_nu_1, ...
                 "Name", satNames{1});

for p = 1:nPlanes
    for s = 1:nSatsPerPlane
        idx = (p-1)*nSatsPerPlane + s;
        if idx == 1, continue; end          % Already created above
        final_raan = raan_base + plane_raan_offsets(p);   % This plane's RAAN
        final_nu   = nu_start_exact + train_nu_offsets(s); % Position in train
        sats(idx) = satellite(sc, a, ecc, base_inc, final_raan, argp, final_nu, ...
                              "Name", satNames{idx});
    end
end
fprintf('    Created %d satellites.\n', numel(sats));

% ---- Record which plane each satellite belongs to ----
% Used later for colour-coding plots (same colour per plane).
satPlaneIdx = zeros(1, nSats);
for p = 1:nPlanes
    for s = 1:nSatsPerPlane
        satPlaneIdx((p-1)*nSatsPerPlane + s) = p;
    end
end

% ---- Create Ground Station objects ----
% These represent the receiver locations (6 cities in Lombardy).
% The Satellite Communications Toolbox can compute access (visibility),
% look-angles, and range between any satellite and ground station.
cityNetwork = {
    'Milan',    45.4642, 9.1900;     % Primary target
    'Monza',    45.5845, 9.2744;     % ~15 km NE of Milan
    'Pavia',    45.1847, 9.1582;     % ~35 km S of Milan
    'Novara',   45.4469, 8.6212;     % ~50 km W of Milan
    'Bergamo',  45.6983, 9.6773;     % ~50 km NE of Milan
    'Lodi',     45.3165, 9.5032;     % ~30 km SE of Milan
};
numCities = size(cityNetwork, 1);    % 6

gs = groundStation(sc, ...
    'Name',      cityNetwork(:,1)', ...
    'Latitude',  [cityNetwork{:,2}], ...
    'Longitude', [cityNetwork{:,3}]);

% Compute access (line-of-sight visibility) intervals between every
% satellite–city pair.  An "access" exists whenever the satellite is
% above the local horizon (elevation > 0°).
for sIdx = 1:nSats
    for cIdx = 1:numCities
        ac(sIdx, cIdx) = access(sats(sIdx), gs(cIdx));
    end
end

% ----- export results into the shared context ------------------------
ctx.sc            = sc;
ctx.sats          = sats;
ctx.gs            = gs;
ctx.ac            = ac;
ctx.satNames      = satNames;
ctx.satPlaneIdx   = satPlaneIdx;
ctx.nSats         = nSats;
ctx.nPlanes       = nPlanes;
ctx.nSatsPerPlane = nSatsPerPlane;
ctx.numCities     = numCities;
ctx.cityNetwork   = cityNetwork;
ctx.wgs84         = wgs84;
ctx.zenithTime    = zenithTime;
ctx.raan_base     = raan_base;
ctx.nu_start_exact= nu_start_exact;

end
