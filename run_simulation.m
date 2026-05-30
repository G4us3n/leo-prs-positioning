%RUN_SIMULATION  Launcher script for the LEO-PRS positioning simulation.
%
%   Run this file from the repository root. It puts the source folders on
%   the MATLAB path and then executes the full pipeline with the default
%   configuration. The final context struct is left in the base workspace
%   as `ctx` for inspection.
%
%   To customise a run, copy the three lines below into your own script:
%       cfg = config();
%       cfg.nPlanes = 6;        % example: change the geometry
%       ctx = main(cfg);

clear; clc; close all force;

thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(thisDir, 'src')));

ctx = main();
