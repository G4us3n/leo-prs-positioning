%RUN_SIMULATION  Launcher script for the LEO-PRS positioning simulation.
%
%   Open LeoPrsPositioning.prj first, then run this file. The project puts
%   the source folders (src, src/stages) on the MATLAB path automatically,
%   so no addpath is needed here. The script executes the full pipeline
%   with the default configuration and leaves the final context struct in
%   the base workspace as `ctx` for inspection.
%
%   To customise a run, copy the three lines below into your own script:
%       cfg = config();
%       cfg.nPlanes = 6;        % example: change the geometry
%       ctx = main(cfg);
%
%   Running without the project open? Call startup() once first to set the
%   path manually.

clear; clc; close all force;

ctx = main();
