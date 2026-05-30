function ctx = stage07_prs_waveform(cfg, ctx)
%STAGE07_PRS_WAVEFORM  Generate the 3GPP 5G-NR PRS waveform (TS 38.211).
%
%   Pipeline stage. Reads CTX/CFG, runs one section, returns updated CTX.

c=ctx.c;

% ----- original section body (unchanged physics) ---------------------
disp('======================================================');
disp('Section 6: 5G PRS Generation (SCS=30 kHz, 20 MHz FR1)');
disp('======================================================');

% ---- Carrier configuration ----
carrier                   = nrCarrierConfig;
carrier.SubcarrierSpacing = cfg.prs.subcarrierSpacing;          % 30 kHz SCS
carrier.NSizeGrid         = cfg.prs.nSizeGrid;          % 51 RBs → ≈ 20 MHz
carrier.CyclicPrefix      = 'normal';    % Normal cyclic prefix
carrier.NSlot             = 0;           % Slot index
carrier.NCellID           = cfg.prs.nCellID;           % Cell ID (affects scrambling)

% ---- PRS configuration ----
prsCfg                          = nrPRSConfig;
prsCfg.PRSResourceSetPeriod     = [10 0];      % Period and offset in slots
prsCfg.PRSResourceOffset        = 0;
prsCfg.PRSResourceRepetition    = 1;
prsCfg.NumRB                    = cfg.prs.nSizeGrid;           % All 51 RBs carry PRS
prsCfg.RBOffset                 = 0;
prsCfg.CombSize                 = cfg.prs.combSize;            % Every 2nd subcarrier
prsCfg.NumPRSSymbols            = cfg.prs.numPRSSymbols;           % 12 OFDM symbols (≈ full slot)
prsCfg.SymbolStart              = 0;
prsCfg.REOffset                 = 0;
prsCfg.NPRSID                   = cfg.prs.nprsID;          % PRS sequence ID

% ---- Generate the waveform ----
prsSym = nrPRS(carrier, prsCfg);            % QPSK symbols
prsInd = nrPRSIndices(carrier, prsCfg);     % Indices into the resource grid
txGrid = nrResourceGrid(carrier);           % Empty time-frequency grid
txGrid(prsInd) = prsSym;                    % Place PRS symbols

% OFDM modulate: freq domain → time domain
[prs_tx, ofdmInfo] = nrOFDMModulate(carrier, txGrid);
prs_tx = prs_tx / sqrt(mean(abs(prs_tx).^2));   % Normalise to unit power

% ---- Key parameters ----
fs_prs           = ofdmInfo.SampleRate;          % Sample rate [Hz]
Ts_prs           = 1 / fs_prs;                   % Sample period [s]
range_per_sample = c * Ts_prs;                    % ≈ 9.76 m
L_prs            = length(prs_tx);                % Waveform length [samples]

SCS_Hz      = carrier.SubcarrierSpacing * 1e3;    % 30 000 Hz
N_active_SC = prsCfg.NumRB * 12;                  % 612 active subcarriers
N_pilot_SC  = N_active_SC / prsCfg.CombSize;      % 306 pilot subcarriers
M_prs       = prsCfg.NumPRSSymbols;               % 12 symbols
BW_prs      = N_active_SC * SCS_Hz;               % 18.36 MHz active bandwidth

fprintf('\n>>> 5G PRS Waveform Parameters <<<\n');
fprintf('  Carrier SCS        : %d kHz\n',       carrier.SubcarrierSpacing);
fprintf('  Channel BW         : 20 MHz  (%d RBs)\n', prsCfg.NumRB);
fprintf('  Active SC (N)      : %d\n',            N_active_SC);
fprintf('  Pilot SC (N/K)     : %d  (CombSize K=%d)\n', N_pilot_SC, prsCfg.CombSize);
fprintf('  PRS Symbols (M)    : %d\n',            M_prs);
fprintf('  Active BW          : %.2f MHz\n',      BW_prs/1e6);
fprintf('  Sample Rate (fs)   : %.4f MHz\n',      fs_prs/1e6);
fprintf('  Range / sample     : %.4f m\n',        range_per_sample);
fprintf('  Waveform length    : %d samples  (%.3f ms)\n', L_prs, L_prs*Ts_prs*1e3);
fprintf('  FFT size (Nfft)    : %d\n',            ofdmInfo.Nfft);

% ----- export results into the shared context ------------------------
ctx.carrier=carrier; ctx.prsCfg=prsCfg; ctx.prs_tx=prs_tx; ctx.ofdmInfo=ofdmInfo;
ctx.fs_prs=fs_prs; ctx.Ts_prs=Ts_prs; ctx.range_per_sample=range_per_sample;
ctx.L_prs=L_prs; ctx.N_active_SC=N_active_SC; ctx.N_pilot_SC=N_pilot_SC;
ctx.M_prs=M_prs; ctx.BW_prs=BW_prs;

end
