function startup()
%STARTUP  Project startup hook for the LEO-PRS positioning simulation.
%
%   When the MATLAB Project is opened, the Project Path entries (src and
%   src/stages) are added automatically, so this file normally has nothing
%   to do. It is kept so the project can also be initialised without the
%   IDE: running it from the project root puts the source folders on the
%   path, mirroring the behaviour the old run_simulation.m used to provide.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'src'));
addpath(fullfile(thisDir, 'src', 'stages'));
end
