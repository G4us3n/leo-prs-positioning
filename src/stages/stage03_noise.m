function ctx = stage03_noise(cfg, ctx)
%STAGE03_NOISE  Establish receiver thermal noise floor and atmospheric loss model.
%
%   Stage of the LEO-PRS positioning pipeline. Reads the shared context
%   struct CTX (and the CFG parameters), runs one logical section of the
%   simulation, and returns CTX with new fields added. See config.m for
%   parameters and main.m for the run order.

BW          = cfg.BW_Hz;
NF          = cfg.NF_dB;
T_ref       = cfg.T_ref_K;
freq_Spare  = ctx.freq_Spare;
centerLat   = cfg.centerLat;
centerLon   = cfg.centerLon;

% ----- original section body (unchanged physics) ---------------------
k_boltzmann = 1.380649e-23;      % Boltzmann constant [J/K]

% kTB noise power  →  convert to dBm  →  add NF
noise_power_W     = k_boltzmann * T_ref * BW;            % ≈ 4.00 × 10⁻¹⁵ W
noise_floor_dBm   = 10*log10(noise_power_W) + 30 + NF;   % ≈ −101 dBm
noise_density_dBm = 10*log10(k_boltzmann * T_ref) + 30 + NF;  % ≈ −171 dBm/Hz

% ---- ITU-R atmospheric-loss model configuration ----
% These parameters feed into p618PropagationLosses() which implements
% the full ITU-R P.618 model: gaseous absorption (O₂, H₂O), rain,
% cloud, and scintillation attenuation along the slant path.
cfg_atm = struct();
cfg_atm.Frequency       = freq_Spare;       % 20 GHz
cfg_atm.ElevationAngle  = 5;                % placeholder (overridden per link)
cfg_atm.Latitude        = centerLat;
cfg_atm.Longitude       = centerLon;
cfg_atm.TotalColumnWaterVapourDensity = cfg.waterVapourDensity;  % kg/m² (typical mid-latitude)
cfg_atm.Polarization    = cfg.polarization;               % Vertical polarisation
cfg_atm.AnnualExceedance = cfg.annualExceedance;                % Design for 1% worst-case rain

fprintf('\n>>> Communication Physical Layer Parameters Established <<<\n');
fprintf('    Bandwidth: 10 MHz\n    System Noise Floor: %.2f dBm\n', noise_floor_dBm);
fprintf('    Noise Power Density: %.2f dBm/Hz\n', noise_density_dBm);
fprintf('    Atmospheric Model: ITU-R P.618 (gas + rain + clouds + scintillation)\n');

% ----- export results into the shared context ------------------------
ctx.BW                = BW;
ctx.NF                = NF;
ctx.T_ref             = T_ref;
ctx.noise_floor_dBm   = noise_floor_dBm;
ctx.noise_density_dBm = noise_density_dBm;
ctx.cfg_atm           = cfg_atm;

end
