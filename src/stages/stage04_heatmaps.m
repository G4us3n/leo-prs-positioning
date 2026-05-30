function ctx = stage04_heatmaps(cfg, ctx)
%STAGE04_HEATMAPS  Geographic coverage heatmaps at two snapshots.
%
%   Pipeline stage. Reads CTX/CFG, runs one section, returns updated CTX.

sats=ctx.sats; gs=ctx.gs; nSats=ctx.nSats; nPlanes=ctx.nPlanes;
nSatsPerPlane=ctx.nSatsPerPlane; numCities=ctx.numCities; cityNetwork=ctx.cityNetwork;
wgs84=ctx.wgs84; zenithTime=ctx.zenithTime; satNames=ctx.satNames;
satPlaneIdx=ctx.satPlaneIdx; freq_Spare=ctx.freq_Spare; c=ctx.c; lambda=ctx.lambda;
pat3D=ctx.pat3D; AzMat=ctx.AzMat; ElMat=ctx.ElMat;
noise_floor_dBm=ctx.noise_floor_dBm; noise_density_dBm=ctx.noise_density_dBm;
T_ref=ctx.T_ref; cfg_atm=ctx.cfg_atm;
centerLat=cfg.centerLat; centerLon=cfg.centerLon;
latlim=cfg.latlim; lonlim=cfg.lonlim; spacing=cfg.spacing;

% ----- original section body (unchanged physics) ---------------------
timeSteps  = [zenithTime - seconds(6), zenithTime + seconds(6)];
timeLabels = {'01-Apr-2026 07:59:54 UTC (0s)', '01-Apr-2026 08:00:06 UTC (+12s)'};

% Geographic grid: Lombardy region, 0.01° resolution ≈ 1.1 km
latlim  = [44.0, 47.0];
lonlim  = [7.4, 11.0];
spacing = 0.01;
latVec  = latlim(1):spacing:latlim(2);
lonVec  = lonlim(1):spacing:lonlim(2);
[LonGrid, LatGrid] = meshgrid(lonVec, latVec);

% Convert entire grid to ECEF (Earth-Centred Earth-Fixed) coordinates.
% ECEF: origin = Earth centre, X → 0° lon, Z → North Pole.
[x_grid, y_grid, z_grid] = geodetic2ecef(wgs84, LatGrid, LonGrid, 0);

% Three metrics to compute
metricConfig = {
%   Name                  Data   cMin   cMax  Step  Unit
    'Combined Rx Power',  [],   -120,  -60,   4,   'dBm';
    'Combined SNR',       [],    -20,   40,   4,   'dB';
    'Combined C/N0',      [],     50,  110,   4,   'dB-Hz'
};

for tIdx = 1:2   % Loop over the two time instants
    currentTime  = timeSteps(tIdx);
    currentLabel = timeLabels{tIdx};

    % Start with −∞ so any real value will "win" the max
    best_rxPowerGrid = -inf(size(LatGrid));
    best_elGrid      = -inf(size(LatGrid));
    satNadirs        = zeros(nSats, 2);         % Sub-satellite points [lat, lon]

    for sIdx = 1:nSats
        % ---- Satellite state (position + velocity) in ECEF ----
        [satPos, satVel] = states(sats(sIdx), currentTime, "CoordinateFrame", "ecef");
        [satLat, satLon, satAlt] = ecef2geodetic(wgs84, satPos(1), satPos(2), satPos(3));
        satNadirs(sIdx, :) = [satLat, satLon];

        % ---- Build the antenna coordinate frame ----
        % The phased array points toward nadir (straight down).
        % We define a local frame:
        %   X_ant → toward nadir (boresight)
        %   Y_ant → cross-track (perpendicular to velocity and nadir)
        %   Z_ant → completes the right-hand system (roughly along-track)
        [nadirX, nadirY, nadirZ] = geodetic2ecef(wgs84, satLat, satLon, 0);
        X_ant = [nadirX; nadirY; nadirZ] - satPos;   % Sat → nadir vector
        X_ant = X_ant / norm(X_ant);                  % Normalise
        Y_ant = cross(X_ant, satVel);                  % Cross-track
        Y_ant = Y_ant / norm(Y_ant);
        Z_ant = cross(X_ant, Y_ant);                  % Along-track

        % ---- Vectors from satellite to every grid point ----
        V_x = x_grid - satPos(1);
        V_y = y_grid - satPos(2);
        V_z = z_grid - satPos(3);
        SlantRange = sqrt(V_x.^2 + V_y.^2 + V_z.^2);  % Distance [m]

        % Unit direction vectors
        V_ux = V_x ./ SlantRange;
        V_uy = V_y ./ SlantRange;
        V_uz = V_z ./ SlantRange;

        % ---- Project into antenna frame → local (az, el) ----
        % Used to look up the antenna gain from pat3D.
        x_proj   = X_ant(1)*V_ux + X_ant(2)*V_uy + X_ant(3)*V_uz;
        y_proj   = Y_ant(1)*V_ux + Y_ant(2)*V_uy + Y_ant(3)*V_uz;
        z_proj   = Z_ant(1)*V_ux + Z_ant(2)*V_uy + Z_ant(3)*V_uz;
        az_local = atan2d(y_proj, x_proj);    % Azimuth in antenna frame
        el_local = asind(z_proj);             % Elevation in antenna frame

        % ---- Ground-station elevation angle ----
        % Angle above the local horizon.  Negative → satellite below horizon.
        R_norm    = sqrt(x_grid.^2 + y_grid.^2 + z_grid.^2);
        sin_el_gs = ((-V_x).*(x_grid./R_norm) + (-V_y).*(y_grid./R_norm) + ...
                     (-V_z).*(z_grid./R_norm)) ./ SlantRange;
        el_gs_deg = asind(sin_el_gs);

        % ---- Atmospheric losses (ITU-R gaseous + rain) ----
        el_gs_clipped = max(el_gs_deg, 5);       % Clip to 5° (model unstable below)

        % gaspl(): gaseous absorption along the slant path.
        % Signature: gaspl(range, freq, T_celsius, P_kPa, rho_wv)
        % It accepts only vectors, so we flatten the grid, compute, reshape.
        T_celsius = T_ref - 273.15;              % 290 K → 16.85 °C
        P_atm_kPa = 101.325;                     % Standard pressure [kPa]
        rho_wv    = 7.5;                          % Water-vapour density [g/m³]
        gridSize  = size(SlantRange);
        atmLoss_vec     = gaspl(SlantRange(:).', freq_Spare, T_celsius, P_atm_kPa, rho_wv);
        atmLossGrid_gas = reshape(atmLoss_vec, gridSize);

        % Rain model: 1.2 dB/km at 20 GHz for moderate rain; rain extends
        % up to 3.5 km (0 °C isotherm at 45°N).
        % Path through rain = height / sin(elevation).
        %rain_specific_dB_km = 1.2;
        %rain_height_km      = 3.5;
        %rain_path_km        = rain_height_km ./ sind(el_gs_clipped);
        %atmLossGrid_rain    = rain_specific_dB_km * rain_path_km;
        atmLossGrid_rain     = 0;

        % Total atmospheric loss (gas + rain) [dB]
        atmLossGrid = atmLossGrid_gas + atmLossGrid_rain;

        % ---- Antenna gain toward each grid point ----
        % interp2 looks up pat3D; points outside the pattern get −60 dBi.
        gainGrid = interp2(AzMat, ElMat, pat3D, az_local, el_local, 'linear', -60);

        % ---- Link budget ----
        prsPower_dBm = 10*log10(cfg.txPower_W) + 30;    % Transmit: 5 W = 36.99 dBm EIRP
        FSPL = reshape(fspl(SlantRange(:).', lambda), gridSize);
            % fspl() computes 20·log₁₀(4πR/λ).  Needs vector input.
        rxGain_dBi = 10*log10(cfg.rxDishEff * (pi * cfg.rxDishDiam_m / lambda)^2);
            % Receive antenna: 0.5 m parabolic dish, 60% efficiency → ≈ 33 dBi

        rxPowerGrid = prsPower_dBm + gainGrid + rxGain_dBi - FSPL - atmLossGrid;

        % ---- Macro-diversity: keep the strongest satellite ----
        best_rxPowerGrid = max(best_rxPowerGrid, rxPowerGrid);
        best_elGrid      = max(best_elGrid, el_gs_deg);
    end

    % Derived metrics
    snrGrid = best_rxPowerGrid - noise_floor_dBm;     % SNR [dB]
    cn0Grid = best_rxPowerGrid - noise_density_dBm;    % C/N₀ [dB-Hz]

    metricConfig{1,2} = best_rxPowerGrid;
    metricConfig{2,2} = snrGrid;
    metricConfig{3,2} = cn0Grid;

    % ---- Plot heatmaps ----
    for mIdx = 1:2    % Rx Power and SNR (skip C/N₀ to save screen space)
        mName = metricConfig{mIdx, 1};
        mData = metricConfig{mIdx, 2};
        mMin  = metricConfig{mIdx, 3};  mMax  = metricConfig{mIdx, 4};
        mStep = metricConfig{mIdx, 5};  mUnit = metricConfig{mIdx, 6};

        figPos = [50 + (mIdx-1)*450, 100 + (tIdx-1)*400, 430, 380];
        figure('Name', sprintf('[%s] %s', currentLabel(1:12), mName), ...
               'Position', figPos, 'Color', 'w');
        worldmap(latlim, lonlim);
        setm(gca, 'FFaceColor', [0.85 0.92 1.0]);        % Light blue = sea
        geoshow('landareas.shp', 'FaceColor', [0.95 0.95 0.90], ...
                'EdgeColor', [0.6 0.6 0.6], 'HandleVisibility', 'off');
        hold on;

        levels = mMin : mStep : mMax;
        colormap(turbo(length(levels)-1));  clim([mMin, mMax]);
        geoshow(LatGrid, LonGrid, mData, "DisplayType","contour", ...
                "Fill","on", "LevelList", levels, "LineColor","none");
        cBar = contourcbar; title(cBar, mUnit, 'FontSize', 10);
        title(sprintf('%s at %s', mName, currentLabel(12:20)), 'FontSize', 11);

        % Annotate each city
        for i = 1:size(cityNetwork, 1)
            cName = cityNetwork{i,1};
            cLat  = cityNetwork{i,2};
            cLon  = cityNetwork{i,3};
            val   = interp2(LonGrid, LatGrid, mData, cLon, cLat, 'linear');
            labelStr = sprintf('%s\n(%.1f)', cName, val);
            if strcmp(cName, 'Milan')
                plotm(cLat, cLon, 'p', 'MarkerSize', 10, ...
                      'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
                textm(cLat-0.015, cLon, labelStr, ...
                      'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', 9);
            else
                plotm(cLat, cLon, 'o', 'MarkerSize', 5, ...
                      'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k');
                textm(cLat-0.015, cLon, labelStr, ...
                      'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', 8);
            end
        end
    end
end

% ----- export results into the shared context ------------------------
% (heatmap stage produces figures only; nothing exported)

end
