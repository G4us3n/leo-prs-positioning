function ctx = stage02_antenna(cfg, ctx)
%STAGE02_ANTENNA  Model the satellite transmit phased array and its radiation pattern.
%
%   Stage of the LEO-PRS positioning pipeline. Reads the shared context
%   struct CTX (and the CFG parameters), runs one logical section of the
%   simulation, and returns CTX with new fields added. See config.m for
%   parameters and main.m for the run order.

freq_Spare = cfg.freq_Hz;
nRows = cfg.arrayRows;  nCols = cfg.arrayCols;

% ----- original section body (unchanged physics) ---------------------
c      = physconst('LightSpeed');           % 299 792 458 m/s
lambda = c / freq_Spare;                    % Wavelength ≈ 15 mm

% Single antenna element: cosine pattern with exponent 1.5.
% This models a patch-type element whose gain falls as cos^1.5(theta)
% from boresight.
elem = phased.CosineAntennaElement('CosinePower', [1.5 1.5]);

% Build the 32 × 32 URA
nRows = 32; nCols = 32;                    % 1024 elements
arrayAnt = phased.URA('Element', elem, ...
    'Size', [nRows nCols], ...
    'ElementSpacing', [0.5 0.5]*lambda, ... % Half-wavelength spacing
    'ArrayNormal', 'x');                    % Array face points along X

% ---- Apply Taylor-window taper + circular mask ----
[X_grid_ant, Y_grid_ant] = meshgrid(1:nCols, 1:nRows);
radiusMap   = sqrt((X_grid_ant - 16.5).^2 + (Y_grid_ant - 16.5).^2);
taperMatrix = taylorwin(nRows, 4, cfg.taperSLL_dB) * taylorwin(nCols, 4, cfg.taperSLL_dB).';
    % taylorwin(N, nbar, SLL):  1-D window with sidelobe level SLL dB.
    % The outer product gives a 2-D taper from two 1-D windows.
taperMatrix(radiusMap > 16.5) = 0;   % Zero out elements outside the circle
arrayAnt.Taper = taperMatrix(:);     % Flatten to column vector

% ---- Compute the full 3-D radiation pattern ----
% Evaluated on a 361 × 361 grid (−90°:0.5°:+90° in both az and el).
% This look-up table is interpolated later to find the gain toward
% any ground point.
az_range = -90:0.5:90;
el_range = -90:0.5:90;
pat3D = pattern(arrayAnt, freq_Spare, az_range, el_range, "PropagationSpeed", c);
[AzMat, ElMat] = meshgrid(az_range, el_range);   % Grids for interp2

% Pre-compute a steering-vector object (for future beam-steering use).
steeringVec = phased.SteeringVector('SensorArray', arrayAnt, ...
    'PropagationSpeed', c);

% ----- export results into the shared context ------------------------
ctx.freq_Spare   = freq_Spare;
ctx.c            = c;
ctx.lambda       = lambda;
ctx.arrayAnt     = arrayAnt;
ctx.pat3D        = pat3D;
ctx.AzMat        = AzMat;
ctx.ElMat        = ElMat;
ctx.steeringVec  = steeringVec;

end
