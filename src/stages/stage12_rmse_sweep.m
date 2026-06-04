function ctx = stage12_rmse_sweep(cfg, ctx)
%STAGE12_RMSE_SWEEP  Ranging RMSE vs SNR sweep for every numerology in cfg.prs_list.
%
%   Runs an independent Monte Carlo sweep for each entry in cfg.prs_list,
%   regenerating the PRS waveform from scratch each time, and overlays all
%   curves on a single figure following the issue-#2 style convention:
%     color  — one per numerology (blue, orange, …)
%     solid  — RMSE with Doppler compensation + parabolic interp
%     dotted — RMSE without Doppler compensation
%     dashed — CRLB (theoretical lower bound)

c = ctx.c;

fprintf('\n>>> Stage 12: RMSE vs SNR — multi-numerology sweep <<<\n');

SNR_sweep_dB    = cfg.mc.snrSweep_dB;
N_MC_sweep      = cfg.mc.nTrialsSweep;
int_delay       = 100;
doppler_test_Hz = cfg.mc.dopplerTest_Hz;
SNR_lin         = 10.^(SNR_sweep_dB / 10);

% ---- Create figure with dual y-axes ----
fig_h = figure('Name', 'Ranging RMSE vs SNR — Multi-Numerology', ...
    'Position', [250 200 1000 650], 'Color', 'w');

yyaxis left;  hold on; grid on; box on;
set(gca, 'YScale', 'log');
ylabel('Range RMSE (m)', 'FontWeight', 'bold', 'FontSize', 11);
ylim([1e-2 1e5]);

yyaxis right; hold on;
ylabel('Detection rate (%)', 'FontWeight', 'bold', 'FontSize', 11);
ylim([0 110]);

all_results = cell(numel(cfg.prs_list), 1);

for n_idx = 1:numel(cfg.prs_list)
    prs_cfg = cfg.prs_list{n_idx};
    col     = prs_cfg.color;

    % ---- Regenerate PRS waveform for this numerology ----
    carrier                   = nrCarrierConfig;
    carrier.SubcarrierSpacing = prs_cfg.subcarrierSpacing;
    carrier.NSizeGrid         = prs_cfg.nSizeGrid;
    carrier.CyclicPrefix      = 'normal';
    carrier.NSlot             = 0;
    carrier.NCellID           = prs_cfg.nCellID;

    prsCfg_n                       = nrPRSConfig;
    prsCfg_n.PRSResourceSetPeriod  = [10 0];
    prsCfg_n.PRSResourceOffset     = 0;
    prsCfg_n.PRSResourceRepetition = 1;
    prsCfg_n.NumRB                 = prs_cfg.nSizeGrid;
    prsCfg_n.RBOffset              = 0;
    prsCfg_n.CombSize              = prs_cfg.combSize;
    prsCfg_n.NumPRSSymbols         = prs_cfg.numPRSSymbols;
    prsCfg_n.SymbolStart           = 0;
    prsCfg_n.REOffset              = 0;
    prsCfg_n.NPRSID                = prs_cfg.nprsID;

    prsSym = nrPRS(carrier, prsCfg_n);
    prsInd = nrPRSIndices(carrier, prsCfg_n);
    txGrid = nrResourceGrid(carrier);
    txGrid(prsInd) = prsSym;
    [prs_tx_n, ofdmInfo] = nrOFDMModulate(carrier, txGrid);
    prs_tx_n = prs_tx_n / sqrt(mean(abs(prs_tx_n).^2));

    Ts_n            = 1 / ofdmInfo.SampleRate;
    L_n             = length(prs_tx_n);
    BW_n            = prs_cfg.nSizeGrid * 12 * prs_cfg.subcarrierSpacing * 1e3;
    range_per_smp_n = c * Ts_n;

    fprintf('\n  [%d] %s | SCS=%d kHz | BW=%.1f MHz | L=%d | range/smp=%.3f m\n', ...
        n_idx, prs_cfg.label, prs_cfg.subcarrierSpacing, BW_n/1e6, L_n, range_per_smp_n);

    % ---- Theoretical bounds ----
    CRLB_n  = c * sqrt(3) ./ (pi * BW_n * sqrt(2 * L_n * SNR_lin));
    quant_n = range_per_smp_n / sqrt(12);

    [snr_fft_dB_n, ~, crlb_fft_chk_n, ~] = snr_after_fft(SNR_sweep_dB, L_n, BW_n, c);
    int_gain_dB_n = 10*log10(L_n);

    % ---- Monte Carlo sweep ----
    RMSE_dop   = nan(numel(SNR_sweep_dB), 1);
    RMSE_nodop = nan(numel(SNR_sweep_dB), 1);
    det_rate   = nan(numel(SNR_sweep_dB), 1);

    for sIdx = 1:numel(SNR_sweep_dB)
        SNRdB    = SNR_sweep_dB(sIdx);
        err_dop  = nan(N_MC_sweep, 1);
        err_nodop = nan(N_MC_sweep, 1);
        detected = 0;

        for mc = 1:N_MC_sweep
            % Random fractional delay each trial (avoids grid-alignment bias)
            frac          = rand();
            ref_delay_smp = int_delay + frac;
            intRef        = floor(ref_delay_smp);
            fracRef       = ref_delay_smp - intRef;

            L_mc     = L_n + intRef + 200;
            Smc      = fft([prs_tx_n; zeros(L_mc - L_n, 1)]);
            kmc      = (0:L_mc-1).';
            kmc(kmc > L_mc/2) = kmc(kmc > L_mc/2) - L_mc;
            sig_frac = ifft(Smc .* exp(-1j*2*pi*kmc*fracRef/L_mc));
            rx_t     = circshift(sig_frac, intRef);
            rx_t(1:intRef) = 0;

            t_vec_mc = (0:L_mc-1).' * Ts_n;
            rx_t_dop = rx_t .* exp(1j*2*pi*doppler_test_Hz*t_vec_mc);

            sig_pow = mean(abs(sig_frac).^2);
            n_pow   = sig_pow / 10^(SNRdB/10);
            noise_v = sqrt(n_pow/2) * (randn(L_mc,1) + 1j*randn(L_mc,1));
            rx_mc   = rx_t_dop + noise_v;

            % With Doppler compensation
            rx_comp = rx_mc .* exp(-1j*2*pi*doppler_test_Hz*t_vec_mc);
            [xcm, lgm] = xcorr(rx_comp, prs_tx_n);
            xcm_pos = abs(xcm(lgm >= 0));  lg_pos = lgm(lgm >= 0);
            [pkv, pki] = max(xcm_pos);
            if pkv >= 4 * median(xcm_pos)
                detected = detected + 1;
                if pki > 1 && pki < numel(xcm_pos)
                    yA = xcm_pos(pki-1); yB = xcm_pos(pki); yC = xcm_pos(pki+1);
                    d  = 0.5*(yA-yC)/(yA-2*yB+yC);
                    pk_l = lg_pos(pki) + d;
                else
                    pk_l = lg_pos(pki);
                end
                err_dop(mc) = (pk_l - ref_delay_smp) * range_per_smp_n;
            end

            % Without Doppler compensation
            [xcm2, lgm2] = xcorr(rx_mc, prs_tx_n);
            xcm2_pos = abs(xcm2(lgm2 >= 0));  lg2_pos = lgm2(lgm2 >= 0);
            [pkv2, pki2] = max(xcm2_pos);
            if pkv2 >= 4 * median(xcm2_pos)
                if pki2 > 1 && pki2 < numel(xcm2_pos)
                    yA2=xcm2_pos(pki2-1); yB2=xcm2_pos(pki2); yC2=xcm2_pos(pki2+1);
                    d2  = 0.5*(yA2-yC2)/(yA2-2*yB2+yC2);
                    pk2 = lg2_pos(pki2) + d2;
                else
                    pk2 = lg2_pos(pki2);
                end
                err_nodop(mc) = (pk2 - ref_delay_smp) * range_per_smp_n;
            end
        end

        valid_dop   = err_dop(~isnan(err_dop));
        valid_nodop = err_nodop(~isnan(err_nodop));
        if ~isempty(valid_dop),   RMSE_dop(sIdx)   = sqrt(mean(valid_dop.^2));   end
        if ~isempty(valid_nodop), RMSE_nodop(sIdx)  = sqrt(mean(valid_nodop.^2)); end
        det_rate(sIdx) = detected / N_MC_sweep;
    end

    % ---- Plot this numerology ----
    col_light = 0.55*col + 0.45*[1 1 1];  % lighter shade for no-Doppler and detection

    yyaxis left;
    semilogy(SNR_sweep_dB, RMSE_dop,   '-',  'Color', col, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('RMSE Doppler comp — %s', prs_cfg.label));
    semilogy(SNR_sweep_dB, RMSE_nodop, ':',  'Color', col_light, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('RMSE no Doppler — %s', prs_cfg.label));
    semilogy(SNR_sweep_dB, CRLB_n,     '--', 'Color', col, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('CRLB — %s', prs_cfg.label));
    if n_idx == 1
        yline(quant_n, 'm:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        text(SNR_sweep_dB(end), quant_n*1.6, 'Sample res.', ...
            'FontSize', 8, 'Color', 'm', 'HorizontalAlignment', 'right');
    end

    yyaxis right;
    plot(SNR_sweep_dB, det_rate*100, ':', 'Color', col_light, 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');

    % ---- Store results ----
    all_results{n_idx} = struct( ...
        'label', prs_cfg.label, 'col', col, ...
        'RMSE_dop', RMSE_dop, 'RMSE_nodop', RMSE_nodop, 'CRLB', CRLB_n, ...
        'det_rate', det_rate, 'L', L_n, 'BW', BW_n, 'quant', quant_n, ...
        'snr_fft_dB', snr_fft_dB_n, 'int_gain_dB', int_gain_dB_n, ...
        'crlb_fft_chk', crlb_fft_chk_n);
end

% ---- Finalise plot ----
yyaxis left;
ax_main = gca;
xlabel('Per-sample SNR  (dB)', 'FontWeight', 'bold', 'FontSize', 11);
title('PRS Ranging Accuracy vs SNR — Multi-Numerology (with Doppler)', ...
    'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 8);

yyaxis right;
ax_main.YColor = [0.4 0.4 0.4];

% ---- Print summary ----
fprintf('\n  %-20s  %-8s  %-8s  %-14s  %-14s\n', ...
    'Numerology', 'BW(MHz)', 'L_prs', 'CRLB@0dB(m)', 'CRLB@20dB(m)');
fprintf('  %s\n', repmat('-', 1, 68));
for n_idx = 1:numel(cfg.prs_list)
    r = all_results{n_idx};
    fprintf('  %-20s  %-8.1f  %-8d  %-14.4f  %-14.4f\n', ...
        r.label, r.BW/1e6, r.L, ...
        r.CRLB(SNR_sweep_dB == 0), r.CRLB(SNR_sweep_dB == 20));
end

fprintf('\n  --- SNR-after-FFT verification (paper §II-B, eq. 7) ---\n');
for n_idx = 1:numel(cfg.prs_list)
    r = all_results{n_idx};
    fprintf('  %s — coherent integration gain: +%.1f dB  (L=%d)\n', ...
        r.label, r.int_gain_dB, r.L);
    fprintf('  %-16s  %-18s  %-12s  %s\n', ...
        'SNR_sample(dB)', 'SNR_post-FFT(dB)', 'CRLB_fft(m)', '|Delta CRLB|(m)');
    fprintf('  %s\n', repmat('-', 1, 68));
    for snr_key = [-10, 0, 10, 20]
        idx = find(SNR_sweep_dB == snr_key, 1);
        if ~isempty(idx)
            diff_m = abs(r.crlb_fft_chk(idx) - r.CRLB(idx));
            fprintf('  %+12.1f dB      %+14.1f dB      %8.4f m    %.2e m\n', ...
                SNR_sweep_dB(idx), r.snr_fft_dB(idx), r.crlb_fft_chk(idx), diff_m);
        end
    end
end

% ---- Export — first numerology kept for backward compat ----
r1 = all_results{1};
ctx.SNR_sweep_dB         = SNR_sweep_dB;
ctx.CRLB_dist            = r1.CRLB;
ctx.range_RMSE_interp    = r1.RMSE_dop;
ctx.range_RMSE_raw       = r1.RMSE_dop;
ctx.range_RMSE_noDopComp = r1.RMSE_nodop;
ctx.detect_rate          = r1.det_rate;
ctx.all_rmse_results     = all_results;

end
