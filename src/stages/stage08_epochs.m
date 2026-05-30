function ctx = stage08_epochs(cfg, ctx)
%STAGE08_EPOCHS  Select reference epochs and build per-epoch best-satellite structs.
%
%   Pipeline stage. Reads CTX/CFG, runs one section, returns updated CTX.

sats=ctx.sats; gs=ctx.gs; nSats=ctx.nSats; nPlanes=ctx.nPlanes;
nSatsPerPlane=ctx.nSatsPerPlane; numCities=ctx.numCities; cityNetwork=ctx.cityNetwork;
wgs84=ctx.wgs84; zenithTime=ctx.zenithTime; satNames=ctx.satNames;
satPlaneIdx=ctx.satPlaneIdx; freq_Spare=ctx.freq_Spare; c=ctx.c; lambda=ctx.lambda;
pat3D=ctx.pat3D; AzMat=ctx.AzMat; ElMat=ctx.ElMat;
noise_floor_dBm=ctx.noise_floor_dBm; noise_density_dBm=ctx.noise_density_dBm;
T_ref=ctx.T_ref; cfg_atm=ctx.cfg_atm;
relTimeMin=ctx.relTimeMin;
history_all_Elev=ctx.history_all_Elev; history_all_Range=ctx.history_all_Range;
history_all_SNR=ctx.history_all_SNR;

% ----- original section body (unchanged physics) ---------------------
milanIdx = 1;   % Index of Milan in cityNetwork

% ---- Find the zenith epoch for each plane ----
planeZenithOffsets = zeros(nPlanes, 1);
planeZenithIdx     = zeros(nPlanes, 1);
planeZenithElev    = zeros(nPlanes, 1);
planeZenithSat     = zeros(nPlanes, 1);

for p = 1:nPlanes
    bestElev = -Inf;  bestTimeIdx = 1;  bestSat = 1;
    for s = 1:nSatsPerPlane
        sIdx = (p-1)*nSatsPerPlane + s;
        elev_track = history_all_Elev(:, milanIdx, sIdx);
        [peakEl, peakIdx] = max(elev_track);
        if ~isnan(peakEl) && peakEl > bestElev
            bestElev = peakEl;  bestTimeIdx = peakIdx;  bestSat = sIdx;
        end
    end
    planeZenithElev(p)    = bestElev;
    planeZenithIdx(p)     = bestTimeIdx;
    planeZenithOffsets(p) = relTimeMin(bestTimeIdx);
    planeZenithSat(p)     = bestSat;
end

% Fixed epochs at ±1.5 min
fixedOffsets = [-1.5, +1.5];
fixedIdx = zeros(numel(fixedOffsets), 1);
for k = 1:numel(fixedOffsets)
    [~, fixedIdx(k)] = min(abs(relTimeMin - fixedOffsets(k)));
end

% Combine, sort, and de-duplicate
allOffsets = [planeZenithOffsets; fixedOffsets(:)];
allIdx     = [planeZenithIdx;    fixedIdx];
allTags    = cell(nPlanes + numel(fixedOffsets), 1);
for p = 1:nPlanes
    allTags{p} = sprintf('Plane-%d Zenith', p);
end
allTags{nPlanes+1} = 'Pre (-1.5 min)';
allTags{nPlanes+2} = 'Post (+1.5 min)';

[allOffsets, sortOrd] = sort(allOffsets);
allIdx  = allIdx(sortOrd);
allTags = allTags(sortOrd);
[~, uniOrd] = unique(allIdx, 'stable');
allOffsets = allOffsets(uniOrd);
allIdx     = allIdx(uniOrd);
allTags    = allTags(uniOrd);

nRef       = numel(allOffsets);
refTimeIdx = allIdx;

% Print summary table
ref_labels = cell(nRef, 1);
for r = 1:nRef
    abs_time = zenithTime + minutes(allOffsets(r));
    ref_labels{r} = sprintf('%s | t = %+.2f min | %s UTC', ...
        allTags{r}, allOffsets(r), datestr(abs_time, 'HH:MM:SS'));
end

fprintf('\n>>> Reference Epochs (%d total, sorted by time) <<<\n', nRef);
fprintf('---------------------------------------------------------------------------\n');
fprintf('%-3s | %-18s | %-10s | %-14s | %-10s\n', '#', 'Tag', 't_rel(min)', 't_abs(UTC)', 'Milan-Elev');
fprintf('---------------------------------------------------------------------------\n');
for r = 1:nRef
    abs_time = zenithTime + minutes(allOffsets(r));
    elev_milan_here = max(history_all_Elev(refTimeIdx(r), milanIdx, :), [], 'all');
    fprintf('%-3d | %-18s | %+8.2f   | %-14s | %6.2f°\n', ...
        r, allTags{r}, allOffsets(r), datestr(abs_time, 'HH:MM:SS'), elev_milan_here);
end

% ---- Build epoch structs ----
ep = repmat(struct('label','','idx',[],'tag','', ...
                   'bestSat',[],'SNR',[],'Range',[],'Elev',[], ...
                   'visMask',[],'RMSE',[],'detectRate',[]), nRef, 1);

for r = 1:nRef
    rIdx    = refTimeIdx(r);
    snr_all = squeeze(history_all_SNR  (rIdx, :, :));      % numCities × nSats
    rng_all = squeeze(history_all_Range(rIdx, :, :)) * 1e3; % km → m
    elv_all = squeeze(history_all_Elev (rIdx, :, :));

    [bestSNR, bestSat] = max(snr_all, [], 2);

    bestRange = nan(numCities, 1);
    bestElev  = nan(numCities, 1);
    for cIdx = 1:numCities
        if ~isnan(bestSNR(cIdx))
            bestRange(cIdx) = rng_all(cIdx, bestSat(cIdx));
            bestElev(cIdx)  = elv_all(cIdx, bestSat(cIdx));
        end
    end

    ep(r).label      = ref_labels{r};
    ep(r).idx        = rIdx;
    ep(r).tag        = allTags{r};
    ep(r).bestSat    = bestSat;
    ep(r).SNR        = bestSNR;
    ep(r).Range      = bestRange;
    ep(r).Elev       = bestElev;
    ep(r).visMask    = ~isnan(bestSNR);
    ep(r).RMSE       = nan(numCities, 1);
    ep(r).detectRate = nan(numCities, 1);
end

% ----- export results into the shared context ------------------------
ctx.ep=ep; ctx.nRef=nRef; ctx.refTimeIdx=refTimeIdx; ctx.allTags=allTags;
ctx.ref_labels=ref_labels; ctx.milanIdx=milanIdx;

end
