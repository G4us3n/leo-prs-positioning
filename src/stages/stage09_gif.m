function ctx = stage09_gif(cfg, ctx)
%STAGE09_GIF  Optional animated GIF of the SNR coverage footprint during the flyover.
%
%   Pipeline stage. Reads CTX/CFG, runs one section, returns updated CTX.

sats=ctx.sats; nSats=ctx.nSats; numCities=ctx.numCities; cityNetwork=ctx.cityNetwork;
wgs84=ctx.wgs84; zenithTime=ctx.zenithTime; satNames=ctx.satNames;
freq_Spare=ctx.freq_Spare; c=ctx.c; lambda=ctx.lambda;
pat3D=ctx.pat3D; AzMat=ctx.AzMat; ElMat=ctx.ElMat;
noise_floor_dBm=ctx.noise_floor_dBm; T_ref=ctx.T_ref; cfg_atm=ctx.cfg_atm;
timeVec=ctx.timeVec; relTimeMin=ctx.relTimeMin;
history_all_Elev=ctx.history_all_Elev; history_all_Range=ctx.history_all_Range;
history_all_SNR=ctx.history_all_SNR; history_all_Doppler=ctx.history_all_Doppler;
prs_tx=ctx.prs_tx; Ts_prs=ctx.Ts_prs; range_per_sample=ctx.range_per_sample;
L_prs=ctx.L_prs; BW_prs=ctx.BW_prs;
ep=ctx.ep; nRef=ctx.nRef; refTimeIdx=ctx.refTimeIdx; allTags=ctx.allTags;
numEpochs=ctx.numEpochs; sat_Trace_Lat=ctx.sat_Trace_Lat;
sat_Trace_Lon=ctx.sat_Trace_Lon; latlim=cfg.latlim; lonlim=cfg.lonlim; spacing=cfg.spacing;
centerLat=cfg.centerLat; centerLon=cfg.centerLon;

% ----- original section body (unchanged physics) ---------------------
fprintf('\n>>> Section 7c: Generating SNR Heatmap Animation (GIF) <<<\n');

% ---- Animation parameters ----
gif_step_s  = 10;       % One frame every 10 seconds (37 frames total)
                         %   Use 5 for smoother animation (73 frames)
                         %   Use 20 for faster computation  (19 frames)
gif_delay_s = 0.3;      % Delay between frames in the GIF [s]
                         %   0.3 gives a comfortable viewing speed
gif_dpi     = 150;      % Resolution of each frame
                         %   150 is a good balance between quality and file size

% ---- Build the time vector for the animation frames ----
gif_timeVec = (zenithTime - seconds(90)) : seconds(gif_step_s) : (zenithTime + seconds(90));
nFrames     = length(gif_timeVec);
fprintf('  Frames to generate: %d  (one every %d s)\n', nFrames, gif_step_s);

% ---- Create output folder for individual frames ----
gifDir = 'heatmap_frames';
if ~exist(gifDir, 'dir')
    mkdir(gifDir);
    fprintf('  Created folder: %s/\n', gifDir);
end

% ---- SNR heatmap colour scale (fixed across all frames for consistency) ----
snr_levels = -20:4:40;    % Contour levels [dB]
snr_cmin   = -20;          % Colour axis minimum [dB]
snr_cmax   = 40;           % Colour axis maximum [dB]

% ---- Generate each frame ----
for fr = 1:nFrames
    currentTime = gif_timeVec(fr);
    relSec      = seconds(currentTime - zenithTime);  % Time relative to zenith [s]

    % --- Compute best-of-all-satellites received power on the geographic grid ---
    %     This is the same link-budget computation used in Section 4 and 7b,
    %     repeated here for each animation frame.
    best_rxPowerGrid = -inf(size(LatGrid));   % Start with -inf so any real value wins

    for sIdx = 1:nSats
        % Satellite state at this frame's time instant
        [satPos, satVel] = states(sats(sIdx), currentTime, "CoordinateFrame", "ecef");
        [satLat, satLon, ~] = ecef2geodetic(wgs84, satPos(1), satPos(2), satPos(3));

        % Nadir point in ECEF (used to build the antenna coordinate frame)
        [nadirX, nadirY, nadirZ] = geodetic2ecef(wgs84, satLat, satLon, 0);

        % Antenna coordinate frame: boresight toward nadir
        X_ant = [nadirX; nadirY; nadirZ] - satPos;  X_ant = X_ant / norm(X_ant);
        Y_ant = cross(X_ant, satVel);                Y_ant = Y_ant / norm(Y_ant);
        Z_ant = cross(X_ant, Y_ant);

        % Vectors from satellite to every grid point
        V_x = x_grid - satPos(1);
        V_y = y_grid - satPos(2);
        V_z = z_grid - satPos(3);
        SlantRange = sqrt(V_x.^2 + V_y.^2 + V_z.^2);

        % Unit direction vectors
        V_ux = V_x ./ SlantRange;
        V_uy = V_y ./ SlantRange;
        V_uz = V_z ./ SlantRange;

        % Project into the antenna frame to get local (azimuth, elevation)
        x_proj   = X_ant(1)*V_ux + X_ant(2)*V_uy + X_ant(3)*V_uz;
        y_proj   = Y_ant(1)*V_ux + Y_ant(2)*V_uy + Y_ant(3)*V_uz;
        z_proj   = Z_ant(1)*V_ux + Z_ant(2)*V_uy + Z_ant(3)*V_uz;
        az_local = atan2d(y_proj, x_proj);
        el_local = asind(z_proj);

        % Ground-station elevation angle (clipped to 5 deg for atmospheric model)
        R_norm    = sqrt(x_grid.^2 + y_grid.^2 + z_grid.^2);
        sin_el_gs = ((-V_x).*(x_grid./R_norm) + (-V_y).*(y_grid./R_norm) + ...
                     (-V_z).*(z_grid./R_norm)) ./ SlantRange;
        el_gs_clipped = max(asind(sin_el_gs), 5);

        % Atmospheric losses: gaseous absorption + rain attenuation
        gridSz    = size(SlantRange);
        atmLoss_v = gaspl(SlantRange(:).', freq_Spare, T_ref - 273.15, 101.325, 7.5);
        atmLossGrid = reshape(atmLoss_v, gridSz) + ...
                      1.2 * (3.5 ./ sind(el_gs_clipped));   % Rain: 1.2 dB/km, 3.5 km height

        % Antenna gain from pre-computed pattern look-up table
        gainGrid = interp2(AzMat, ElMat, pat3D, az_local, el_local, 'linear', -60);

        % Link budget: Tx power + Tx gain + Rx gain - FSPL - atmospheric loss
        prsPower_dBm = 10*log10(cfg.txPower_W) + 30;                                    % 5 W EIRP
        FSPL         = reshape(fspl(SlantRange(:).', lambda), gridSz);       % Free-space path loss
        rxGain_dBi   = 10*log10(cfg.rxDishEff * (pi * cfg.rxDishDiam_m / lambda)^2);              % 0.5 m dish, 60% eff.

        rxPowerGrid = prsPower_dBm + gainGrid + rxGain_dBi - FSPL - atmLossGrid;

        % Macro-diversity: keep the strongest satellite at each grid point
        best_rxPowerGrid = max(best_rxPowerGrid, rxPowerGrid);
    end

    % Convert received power to SNR
    snrGrid = best_rxPowerGrid - noise_floor_dBm;

% --- Plot the heatmap for this frame ---
    fig = figure('Visible', 'off', 'Position', [100 100 700 550], 'Color', 'w');

    gx = geoaxes(fig);

    % Use satellite basemap for terrain imagery
    try
        geobasemap(gx, 'satellite');
    catch
        try
            geobasemap(gx, 'landcover');
        catch
            geobasemap(gx, 'streets-light');
        end
    end

    hold(gx, 'on');
    geolimits(gx, latlim, lonlim);

    % --- Plot SNR as a dense scatter overlay ---
    %     We subsample the grid to keep rendering fast while still
    %     producing a smooth-looking heatmap.  Every 3rd point in each
    %     direction gives ~40 000 dots, which fills the map nicely.
    step = 3;
    lat_sub = LatGrid(1:step:end, 1:step:end);
    lon_sub = LonGrid(1:step:end, 1:step:end);
    snr_sub = snrGrid(1:step:end, 1:step:end);

    % Flatten to column vectors for geoscatter
    lat_v = lat_sub(:);
    lon_v = lon_sub(:);
    snr_v = snr_sub(:);

    % Remove points with very low SNR (they would cover terrain uselessly)
    valid = snr_v > (snr_cmin + 2);
    lat_v = lat_v(valid);
    lon_v = lon_v(valid);
    snr_v = snr_v(valid);

    % Plot as filled scatter — marker size controls "pixel" density
    geoscatter(gx, lat_v, lon_v, 4, snr_v, 'filled', ...
               'MarkerFaceAlpha', 0.6);

    % Colormap and colour axis
    colormap(gx, turbo(256));
    clim(gx, [snr_cmin snr_cmax]);
    cb = colorbar(gx, 'eastoutside');
    cb.Label.String = 'SNR (dB)';
    cb.Label.FontSize = 11;

    % --- Mark all ground stations ---
    for cIdx = 1:numCities
        cLat  = cityNetwork{cIdx, 2};
        cLon  = cityNetwork{cIdx, 3};
        cName = cityNetwork{cIdx, 1};

        if strcmp(cName, 'Milan')
            geoplot(gx, cLat, cLon, 'p', 'MarkerSize', 14, ...
                    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
        else
            geoplot(gx, cLat, cLon, 'o', 'MarkerSize', 7, ...
                    'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        end

        text(gx, cLat - 0.08, cLon, cName, ...
             'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
             'FontSize', 8, 'Color', 'w', ...
             'BackgroundColor', [0 0 0 0.5], 'Margin', 1);
    end

    % --- Title ---
    title(gx, sprintf('Combined SNR — t = %+.0f s from zenith', relSec), ...
          'FontSize', 13, 'FontWeight', 'bold');

    % --- Save frame ---
    framePath = fullfile(gifDir, sprintf('frame_%03d.png', fr));
    exportgraphics(fig, framePath, 'Resolution', gif_dpi);
    close(fig);
    fprintf('  Frame %3d/%d  |  t = %+4.0f s  |  saved\n', fr, nFrames, relSec);
end

% ---- Assemble all frames into an animated GIF ----
fprintf('\n  Assembling GIF...\n');
gifName = 'fullfile(cfg.resultsDir,'snr_coverage.gif')';

for fr = 1:nFrames
    framePath = fullfile(gifDir, sprintf('frame_%03d.png', fr));
    [imgRGB, ~] = imread(framePath);

    % Convert RGB image to indexed colour (required by GIF format)
    [imInd, colorMap] = rgb2ind(imgRGB, 256);

    if fr == 1
        % First frame: create the GIF file
        imwrite(imInd, colorMap, gifName, 'gif', ...
                'Loopcount', inf, ...        % Loop forever
                'DelayTime', gif_delay_s);
    else
        % Subsequent frames: append to the existing GIF
        imwrite(imInd, colorMap, gifName, 'gif', ...
                'WriteMode', 'append', ...
                'DelayTime', gif_delay_s);
    end
end

fprintf('  GIF saved: %s  (%d frames, %.1f s playback)\n', ...
        gifName, nFrames, nFrames * gif_delay_s);
fprintf('  Individual frames kept in: %s/\n', gifDir);

% ----- export results into the shared context ------------------------
% (GIF stage writes snr_coverage.gif to results/)

end
