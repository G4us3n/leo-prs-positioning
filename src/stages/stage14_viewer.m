function ctx = stage14_viewer(cfg, ctx)
%STAGE14_VIEWER  Optional interactive 3D satellite scenario viewer.
%
%   Pipeline stage. Reads CTX/CFG, runs one section, returns updated CTX.

sc=ctx.sc; gs=ctx.gs; numCities=ctx.numCities;
centerLat=cfg.centerLat; centerLon=cfg.centerLon;

% ----- original section body (unchanged physics) ---------------------
fprintf('\n>>> Section 12: Interactive 3D Scenario Viewer <<<\n');
try
    basemapList = {'streets', 'satellite', 'topographic', 'landwater', 'darkwater'};
    launched = false;
    for bIdx = 1:numel(basemapList)
        try
            v = satelliteScenarioViewer(sc, 'ShowDetails', true, ...
                'Basemap', basemapList{bIdx});
            fprintf('  Basemap: %s\n', basemapList{bIdx});
            launched = true;
            break;
        catch
            continue;
        end
    end
    if ~launched
        v = satelliteScenarioViewer(sc, 'ShowDetails', true);
        fprintf('  Basemap: default (darkwater — no land visible without internet)\n');
        fprintf('  TIP: Enable internet or install a basemap via Add-On Explorer\n');
        fprintf('       to see land masses on the globe.\n');
    end

    for cIdx = 1:numCities
        gs(cIdx).MarkerColor = [1 0.8 0];   % Yellow
        gs(cIdx).MarkerSize  = 8;
    end
    gs(1).MarkerColor = [1 0 0];             % Milan in red
    gs(1).MarkerSize  = 12;

    campos(v, centerLat, centerLon, 2e6);    % Camera over Milan, 2000 km altitude
    disp('  3D Scenario Viewer launched. Close the window when done.');
catch ME
    fprintf('  Viewer not available: %s\n', ME.message);
end

% ----- export results into the shared context ------------------------
% (3D viewer stage; interactive only)

end
