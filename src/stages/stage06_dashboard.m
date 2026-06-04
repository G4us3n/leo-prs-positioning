function ctx = stage06_dashboard(cfg, ctx)
%STAGE06_DASHBOARD  Ground-track map, Milan dashboard and per-city SNR envelopes.
%
%   Pipeline stage. Reads CTX/CFG, runs one section, returns updated CTX.

sats=ctx.sats; gs=ctx.gs; nSats=ctx.nSats; nPlanes=ctx.nPlanes;
nSatsPerPlane=ctx.nSatsPerPlane; numCities=ctx.numCities; cityNetwork=ctx.cityNetwork;
wgs84=ctx.wgs84; zenithTime=ctx.zenithTime; satNames=ctx.satNames;
satPlaneIdx=ctx.satPlaneIdx; freq_Spare=ctx.freq_Spare; c=ctx.c; lambda=ctx.lambda;
pat3D=ctx.pat3D; AzMat=ctx.AzMat; ElMat=ctx.ElMat;
noise_floor_dBm=ctx.noise_floor_dBm; noise_density_dBm=ctx.noise_density_dBm;
T_ref=ctx.T_ref; cfg_atm=ctx.cfg_atm;
relTimeMin=ctx.relTimeMin; numEpochs=ctx.numEpochs;
history_all_Elev=ctx.history_all_Elev; history_all_Range=ctx.history_all_Range;
history_all_Power=ctx.history_all_Power; history_all_SNR=ctx.history_all_SNR;
history_all_Doppler=ctx.history_all_Doppler;
sat_Trace_Lat=ctx.sat_Trace_Lat; sat_Trace_Lon=ctx.sat_Trace_Lon;

% ----- original section body (unchanged physics) ---------------------
planeColors = [1 0 0; 0 0.7 0; 0 0.3 1; 0.8 0.5 0];
trainStyles = {'-', '--', ':'};

% ---------- Ground-track map ----------
figure('Name', '6-Min Orbital Trajectories (12 Sats)', ...
       'Position', [100, 150, 800, 650], 'Color', 'w');
worldmap([25 65], [-15 35]);
setm(gca, 'MapProjection', 'mercator', 'FFaceColor', [0.9 0.93 0.95]);
geoshow('landareas.shp', 'FaceColor', [0.95 0.95 0.90], 'HandleVisibility', 'off');
hold on;

h_planes = gobjects(1, nPlanes);
for p = 1:nPlanes
    for s = 1:nSatsPerPlane
        sIdx = (p-1)*nSatsPerPlane + s;
        hv = 'off'; if s == 1, hv = 'on'; end   % Only 1st sat per plane in legend
        h = plotm(sat_Trace_Lat(:,sIdx), sat_Trace_Lon(:,sIdx), trainStyles{s}, ...
            'Color', planeColors(p,:), 'LineWidth', 1.8, ...
            'DisplayName', sprintf('Plane %d', p), 'HandleVisibility', hv);
        if s == 1, h_planes(p) = h; end
        % Mark start (○) and end (▲) of track
        plotm(sat_Trace_Lat(1,sIdx),   sat_Trace_Lon(1,sIdx),   'o', ...
              'Color', planeColors(p,:), 'MarkerFaceColor', 'w', ...
              'MarkerSize', 4, 'HandleVisibility', 'off');
        plotm(sat_Trace_Lat(end,sIdx), sat_Trace_Lon(end,sIdx), '^', ...
              'Color', planeColors(p,:), 'MarkerFaceColor', planeColors(p,:), ...
              'MarkerSize', 4, 'HandleVisibility', 'off');
    end
end
milanLat = cityNetwork{1,2};  milanLon = cityNetwork{1,3};
h_milan = plotm(milanLat, milanLon, 'p', 'Color', 'k', 'MarkerFaceColor', 'y', ...
    'MarkerSize', 14, 'LineWidth', 1.5, 'DisplayName', 'Milan GS');
legend([h_planes, h_milan], 'Location', 'southwest', 'FontSize', 9);
title(sprintf('LEO Constellation: %d Satellites, %d Planes (6-Minute Window)', nSats, nPlanes));

% ---------- Milan 6-panel dashboard ----------
figure('Name', 'Milan GS Multi-Link Dynamics (12 Sats)', ...
       'Position', [250, 50, 1400, 750], 'Color', 'w');
plotZenithLine = @() xline(0, 'k--', 'Zenith', ...
    'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment', 'bottom');

% Panel 1: Elevation angle
subplot(2,3,1); hold on; grid on;
for sIdx = 1:nSats
    p = satPlaneIdx(sIdx); s = mod(sIdx-1, nSatsPerPlane)+1;
    plot(relTimeMin, squeeze(history_all_Elev(:,1,sIdx)), trainStyles{s}, ...
         'LineWidth', 1.5, 'Color', planeColors(p,:), 'HandleVisibility', 'off');
end
plotZenithLine(); title('Milan: Elevation Angle');
xlabel('Time (min)'); ylabel('Elevation (°)'); ylim([0 90]); xlim([-3 3]);

% Panel 2: Slant range
subplot(2,3,2); hold on; grid on;
for sIdx = 1:nSats
    p = satPlaneIdx(sIdx); s = mod(sIdx-1, nSatsPerPlane)+1;
    plot(relTimeMin, squeeze(history_all_Range(:,1,sIdx)), trainStyles{s}, ...
         'LineWidth', 1.5, 'Color', planeColors(p,:), 'HandleVisibility', 'off');
end
plotZenithLine(); title('Milan: Slant Range (km)');
xlabel('Time (min)'); ylabel('Distance (km)'); xlim([-3 3]);

% Panel 3: Received power
subplot(2,3,3); hold on; grid on;
for sIdx = 1:nSats
    p = satPlaneIdx(sIdx); s = mod(sIdx-1, nSatsPerPlane)+1;
    plot(relTimeMin, squeeze(history_all_Power(:,1,sIdx)), trainStyles{s}, ...
         'LineWidth', 1.5, 'Color', planeColors(p,:), 'HandleVisibility', 'off');
end
plotZenithLine(); title('Milan: Received Power (dBm)');
xlabel('Time (min)'); ylabel('Power (dBm)'); ylim([-200 -50]); xlim([-3 3]);

% Panel 4: SNR with best-of-12 envelope
subplot(2,3,4); hold on; grid on;
for sIdx = 1:nSats
    p = satPlaneIdx(sIdx); s = mod(sIdx-1, nSatsPerPlane)+1;
    plot(relTimeMin, squeeze(history_all_SNR(:,1,sIdx)), trainStyles{s}, ...
         'LineWidth', 1.5, 'Color', planeColors(p,:), 'HandleVisibility', 'off');
end
snr_best_milan = max(squeeze(history_all_SNR(:,1,:)), [], 2);
plot(relTimeMin, snr_best_milan, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Best of 12');
plotZenithLine(); yline(0, 'r-', '0 dB', 'LineWidth', 1.5, 'HandleVisibility','off');
title('Milan: SNR (dB)');
xlabel('Time (min)'); ylabel('SNR (dB)'); ylim([-100 50]); xlim([-3 3]);
legend('Location', 'southwest', 'FontSize', 8);

% Panel 5: Doppler shift
subplot(2,3,5); hold on; grid on;
for sIdx = 1:nSats
    p = satPlaneIdx(sIdx); s = mod(sIdx-1, nSatsPerPlane)+1;
    plot(relTimeMin, squeeze(history_all_Doppler(:,1,sIdx))/1e3, trainStyles{s}, ...
         'LineWidth', 1.2, 'Color', planeColors(p,:), 'HandleVisibility', 'off');
end
plotZenithLine(); yline(0, 'k:', 'HandleVisibility','off');
title('Milan: Doppler Shift (kHz)');
xlabel('Time (min)'); ylabel('Doppler (kHz)'); xlim([-3 3]);

% Panel 6: Visibility timeline
subplot(2,3,6); hold on; grid on;
for sIdx = 1:nSats
    p = satPlaneIdx(sIdx);
    vis_mask = ~isnan(squeeze(history_all_Elev(:,1,sIdx)));
    fill_x = relTimeMin(vis_mask);
    if ~isempty(fill_x)
        bar_y = ones(size(fill_x)) * sIdx;
        plot(fill_x, bar_y, '.', 'Color', planeColors(p,:), ...
             'MarkerSize', 6, 'HandleVisibility', 'off');
    end
end
plotZenithLine();
title('Milan: Satellite Visibility');
xlabel('Time (min)'); ylabel('Satellite #');
ylim([0.5 nSats+0.5]); xlim([-3 3]);

sgtitle(sprintf('Milan GS — %d-Satellite Multi-Link Dashboard', nSats), 'FontWeight','bold');


%% =====================================================================
%  SECTION 5b: PER-CITY SNR ENVELOPE
%  =====================================================================
%  One subplot per city showing the best-of-12 SNR over time.
%  =====================================================================

figure('Name', 'City SNR Time-Series', 'Position', [100 50 1200 850], 'Color', 'w');
nCols_plot = 3;  nRows_plot = ceil(numCities / nCols_plot);

for cIdx = 1:numCities
    subplot(nRows_plot, nCols_plot, cIdx); hold on; grid on;
    snr_best = max(squeeze(history_all_SNR(:,cIdx,:)), [], 2);
    plot(relTimeMin, snr_best, 'LineWidth', 1.8, 'DisplayName', 'Best-sat');
    xline(0, 'k:', 'LineWidth', 1);
    yline(0, 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    title(cityNetwork{cIdx,1}, 'FontWeight', 'bold');
    xlabel('Time (min)'); ylabel('SNR (dB)'); xlim([-3 3]);
    if cIdx == 1, legend('Location', 'southwest', 'FontSize', 8); end
end
sgtitle('Per-City SNR vs Time — Selection Envelope', 'FontWeight','bold', 'FontSize', 12);

% ----- export results into the shared context ------------------------
% (dashboard stage produces figures only)

end
