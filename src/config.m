function cfg = config()
%CONFIG  Central parameter set for the LEO-PRS positioning simulation.
%
%   cfg = CONFIG() returns a struct holding every tunable parameter used by
%   the simulation pipeline. Edit values here rather than inside the stage
%   functions, so that a single source of truth drives the whole run.
%
%   The fields are grouped by stage. Most of the upcoming review tasks
%   (multi-numerology overlays, alternative constellation geometries,
%   multi-satellite least-squares positioning) are controlled from this
%   file alone.

% ---------------------------------------------------------------------
% Scenario time window
% ---------------------------------------------------------------------
% zenithTime is the instant when the constellation converges overhead the
% target city. The simulation spans +/- windowMinutes around it.
cfg.zenithTime    = datetime(2026, 4, 1, 8, 0, 0);   % 1 April 2026, 08:00 UTC
cfg.windowMinutes = 3;                               % +/- 3 min -> 6 min total
cfg.sampleTime    = 1;                               % scenario propagator step [s]
cfg.flyoverStep   = seconds(2);                      % time-series resolution

% ---------------------------------------------------------------------
% Target location (receiver reference) and ground-station network
% ---------------------------------------------------------------------
cfg.centerLat = 45.4642;   % Milan latitude  [deg N]
cfg.centerLon = 9.1900;    % Milan longitude [deg E]

% City network: {name, latitude, longitude}. Row 1 is the primary target.
cfg.cityNetwork = {
    'Milan',    45.4642, 9.1900;
    'Monza',    45.5845, 9.2744;
    'Pavia',    45.1847, 9.1582;
    'Novara',   45.4469, 8.6212;
    'Bergamo',  45.6983, 9.6773;
    'Lodi',     45.3165, 9.5032;
};

% ---------------------------------------------------------------------
% Constellation geometry  (the "train" baseline: 4 planes x 3 sats)
% ---------------------------------------------------------------------
% These four fields are the main knobs for the "test different satellite
% configurations" task. Increase the train_nu spacing to widen along-track
% separation, or widen plane_raan_offsets to extend transversal coverage.
cfg.nPlanes            = 4;
cfg.nSatsPerPlane      = 3;
cfg.altitude_m         = 6871000;          % semi-major axis [m] -> ~493 km alt
cfg.eccentricity       = 0;
cfg.argPerigee         = 0;
cfg.inclination_deg    = 53.0;
cfg.plane_raan_offsets = [-1, 0.0, +1, +2];   % per-plane RAAN offset [deg]
cfg.train_nu_offsets   = [1.5, 0.0, -1.5];    % per-sat true-anomaly offset [deg]

% ---------------------------------------------------------------------
% Transmit antenna (phased array) and RF front end
% ---------------------------------------------------------------------
cfg.freq_Hz        = 20e9;     % Ka-band downlink carrier [Hz]
cfg.arrayRows      = 32;
cfg.arrayCols      = 32;
cfg.taperSLL_dB    = -20;      % Taylor-window sidelobe level
cfg.txPower_W      = 5;        % transmit power [W]
cfg.rxDishDiam_m   = 0.5;      % receive parabolic dish diameter [m]
cfg.rxDishEff      = 0.6;      % receive dish efficiency

% ---------------------------------------------------------------------
% Receiver noise and atmosphere
% ---------------------------------------------------------------------
cfg.BW_Hz   = 10e6;    % receiver bandwidth [Hz]
cfg.NF_dB   = 3;       % receiver noise figure [dB]
cfg.T_ref_K = 290;     % reference temperature [K]
cfg.waterVapourDensity = 25;   % g/m^2, ITU-R P.618 input
cfg.polarization       = 'V';
cfg.annualExceedance   = 1;    % % worst-case rain for the link budget

% ---------------------------------------------------------------------
% Geographic heatmap grid (Lombardy)
% ---------------------------------------------------------------------
cfg.latlim  = [44.0, 47.0];
cfg.lonlim  = [7.4, 11.0];
cfg.spacing = 0.01;            % grid resolution [deg] ~ 1.1 km

% ---------------------------------------------------------------------
% 5G PRS waveform (TS 38.211) - FR1 baseline numerology (mu = 1, 30 kHz)
% ---------------------------------------------------------------------
cfg.prs.subcarrierSpacing = 30;    % kHz (FR1)
cfg.prs.nSizeGrid         = 51;    % resource blocks (~20 MHz)
cfg.prs.combSize          = 2;
cfg.prs.numPRSSymbols     = 12;
cfg.prs.nCellID           = 1;
cfg.prs.nprsID            = 100;

% ---------------------------------------------------------------------
% Monte Carlo ranging
% ---------------------------------------------------------------------
cfg.mc.nTrials        = 250;    % MC iterations per link (Section 8)
cfg.mc.accTarget_m    = 30;     % pass/fail ranging accuracy target [m]
cfg.mc.detectionGate  = 4;      % peak must exceed gate x median to count
cfg.mc.snrSweep_dB    = -50:2:50;  % SNR sweep for the RMSE curve (Section 10)
cfg.mc.nTrialsSweep   = 200;
cfg.mc.dopplerTest_Hz = 300e3;  % representative Doppler for the sweep

% ---------------------------------------------------------------------
% Output / run control
% ---------------------------------------------------------------------
% Toggle expensive or optional stages here. Disabling the GIF and the 3D
% viewer makes a headless / fast run.
cfg.run.heatmaps        = true;
cfg.run.flyover         = true;
cfg.run.dashboards      = true;
cfg.run.prsWaveform     = true;
cfg.run.referenceEpochs = true;
cfg.run.gif             = false;   % expensive; off by default
cfg.run.monteCarlo      = true;
cfg.run.xcorrPlots      = true;
cfg.run.rmseSweep       = true;
cfg.run.dop             = true;
cfg.run.viewer3D        = false;   % needs a display; off by default

cfg.resultsDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results');
end
