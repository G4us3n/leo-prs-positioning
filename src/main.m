function ctx = main(cfg)
%MAIN  Entry point that runs the full LEO-PRS positioning pipeline.
%
%   ctx = MAIN()    runs with the default configuration from config().
%   ctx = MAIN(cfg) runs with a caller-supplied configuration struct.
%
%   The pipeline is a sequence of stages. Each stage is a function
%   ctx = stageNN_name(cfg, ctx) that reads the shared context CTX, runs
%   one logical section of the simulation, and returns CTX with new fields
%   added. Stages can be toggled on/off through cfg.run.*.
%
%   Returns the final context CTX, which holds every intermediate result
%   (satellite handles, flyover histories, PRS waveform, epoch structs,
%   Monte Carlo RMSE, RMSE-vs-SNR sweep, etc.) for further inspection.
%
%   Example:
%       ctx = main();                 % default run
%       cfg = config();
%       cfg.nPlanes = 6;              % try a different geometry
%       cfg.run.gif = true;          % also render the coverage GIF
%       ctx = main(cfg);

if nargin < 1 || isempty(cfg)
    cfg = config();
end

% Make sure the results directory exists.
if ~exist(cfg.resultsDir, 'dir')
    mkdir(cfg.resultsDir);
end

t0 = tic;
banner('LEO-PRS POSITIONING SIMULATION - START');

ctx = struct();

% ---- Always-on stages: build the world and the signal -----------------
ctx = stage01_constellation(cfg, ctx);
ctx = stage02_antenna(cfg, ctx);
ctx = stage03_noise(cfg, ctx);

% ---- Optional / heavy visual stages -----------------------------------
if cfg.run.heatmaps
    ctx = stage04_heatmaps(cfg, ctx);
end

% ---- Time-series flyover is needed by most later stages ---------------
if cfg.run.flyover
    ctx = stage05_flyover(cfg, ctx);
else
    warning('main:flyoverOff', ...
        'Flyover stage is disabled; epoch/MC/DOP stages will be skipped.');
end

if cfg.run.dashboards && isfield(ctx, 'history_all_SNR')
    ctx = stage06_dashboard(cfg, ctx);
end

if cfg.run.prsWaveform
    ctx = stage07_prs_waveform(cfg, ctx);
end

if cfg.run.referenceEpochs && isfield(ctx, 'history_all_SNR')
    ctx = stage08_epochs(cfg, ctx);
end

if cfg.run.gif && isfield(ctx, 'ep')
    ctx = stage09_gif(cfg, ctx);
end

if cfg.run.monteCarlo && isfield(ctx, 'ep') && isfield(ctx, 'prs_tx')
    ctx = stage10_montecarlo(cfg, ctx);
end

if cfg.run.xcorrPlots && isfield(ctx, 'ep') && isfield(ctx, 'prs_tx')
    ctx = stage11_xcorr(cfg, ctx);
end

if cfg.run.rmseSweep && isfield(ctx, 'prs_tx')
    ctx = stage12_rmse_sweep(cfg, ctx);
end

if cfg.run.dop && isfield(ctx, 'ep')
    ctx = stage13_dop(cfg, ctx);
end

if cfg.run.leastSquares && isfield(ctx, 'history_all_Elev')
    ctx = stage15_least_squares(cfg, ctx);
end

if cfg.run.viewer3D
    ctx = stage14_viewer(cfg, ctx);
end

banner(sprintf('SIMULATION COMPLETE in %.1f s', toc(t0)));
end

% =====================================================================
function banner(msg)
line = repmat('=', 1, 54);
fprintf('\n%s\n%s\n%s\n', line, msg, line);
end
