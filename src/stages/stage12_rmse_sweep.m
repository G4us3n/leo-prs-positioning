function ctx = stage12_rmse_sweep(cfg, ctx)
%STAGE12_RMSE_SWEEP  Ranging RMSE vs SNR sweep with CRLB and no-Doppler-compensation curves.
%
%   Pipeline stage. Reads CTX/CFG, runs one section, returns updated CTX.

c=ctx.c; prs_tx=ctx.prs_tx; Ts_prs=ctx.Ts_prs; range_per_sample=ctx.range_per_sample;
L_prs=ctx.L_prs; BW_prs=ctx.BW_prs;

% ----- original section body (unchanged physics) ---------------------
fprintf('\n>>> Section 10: RMSE vs SNR — with Doppler + parabolic interp <<<\n');

SNR_sweep_dB    = cfg.mc.snrSweep_dB;            % SNR points to evaluate
N_MC_sweep      = cfg.mc.nTrialsSweep;                  % MC iterations per SNR point
int_delay       = 100;                  % Fixed integer delay [samples]
doppler_test_Hz = cfg.mc.dopplerTest_Hz;                % 300 kHz (typical LEO at 20 GHz)

% ---- Theoretical bounds ----
SNR_lin       = 10.^(SNR_sweep_dB / 10);        % Linear SNR
CRLB_dist     = c * sqrt(3) ./ (pi * BW_prs * sqrt(2 * L_prs * SNR_lin));
    % ↑ per-sample CRLB: BW_prs = full active bandwidth; L_prs = integration length
quant_floor_m = range_per_sample / sqrt(12);     % Quantisation-limited RMSE

% Pre-allocate results
range_RMSE_interp    = nan(numel(SNR_sweep_dB), 1);   % With Doppler comp + parabolic interp
range_RMSE_raw       = nan(numel(SNR_sweep_dB), 1);   % With Doppler comp, integer sample only
range_RMSE_noDopComp = nan(numel(SNR_sweep_dB), 1);   % WITHOUT Doppler compensation
detect_rate          = nan(numel(SNR_sweep_dB), 1);

for sIdx = 1:numel(SNR_sweep_dB)
    SNRdB      = SNR_sweep_dB(sIdx);
    err_interp = nan(N_MC_sweep, 1);
    err_raw    = nan(N_MC_sweep, 1);
    err_noDop  = nan(N_MC_sweep, 1);
    detected   = 0;

    for mc = 1:N_MC_sweep
        % ---- Random fractional delay each trial ----
        % This prevents the results from being biased by a specific
        % "lucky" or "unlucky" delay alignment with the sample grid.
        frac          = rand();
        ref_delay_smp = int_delay + frac;
        intRef        = floor(ref_delay_smp);
        fracRef       = ref_delay_smp - intRef;

        % ---- Build delayed reference ----
        L_mc     = L_prs + intRef + 200;
        Smc      = fft([prs_tx; zeros(L_mc - L_prs, 1)]);
        kmc      = (0:L_mc-1).';
        kmc(kmc > L_mc/2) = kmc(kmc > L_mc/2) - L_mc;
        sig_frac = ifft(Smc .* exp(-1j*2*pi*kmc*fracRef/L_mc));
        rx_t     = circshift(sig_frac, intRef);
        rx_t(1:intRef) = 0;

        % ---- Apply Doppler shift ----
        t_vec_mc = (0:L_mc-1).' * Ts_prs;
        rx_t_dop = rx_t .* exp(1j*2*pi*doppler_test_Hz*t_vec_mc);

        % ---- Add AWGN (per-sample SNR convention) ----
        % The noise power is set so that:
        %   SNR_per_sample = mean(|signal|²) / mean(|noise|²)
        sig_pow = mean(abs(sig_frac).^2);
        n_pow   = sig_pow / 10^(SNRdB/10);
        noise_v = sqrt(n_pow/2) * (randn(L_mc,1) + 1j*randn(L_mc,1));
        rx_mc   = rx_t_dop + noise_v;

        % ---- WITH Doppler compensation ----
        rx_mc_comp = rx_mc .* exp(-1j*2*pi*doppler_test_Hz*t_vec_mc);

        [xcm, lgm] = xcorr(rx_mc_comp, prs_tx);
        xcm_pos = abs(xcm(lgm >= 0));
        lg_pos  = lgm(lgm >= 0);
        [pkv, pki] = max(xcm_pos);

        if pkv >= 4 * median(xcm_pos)
            detected = detected + 1;

            % Raw (integer-sample) estimate
            pk_raw      = lg_pos(pki);
            err_raw(mc) = (pk_raw - ref_delay_smp) * range_per_sample;

            % Parabolic interpolation
            if pki > 1 && pki < numel(xcm_pos)
                yA = xcm_pos(pki-1);  yB = xcm_pos(pki);  yC = xcm_pos(pki+1);
                d  = 0.5 * (yA - yC) / (yA - 2*yB + yC);
                pk_interp = lg_pos(pki) + d;
            else
                pk_interp = lg_pos(pki);
            end
            err_interp(mc) = (pk_interp - ref_delay_smp) * range_per_sample;
        end

        % ---- WITHOUT Doppler compensation (to quantify degradation) ----
        [xcm2, lgm2] = xcorr(rx_mc, prs_tx);     % Correlate WITHOUT Doppler removal
        xcm2_pos = abs(xcm2(lgm2 >= 0));
        lg2_pos  = lgm2(lgm2 >= 0);
        [pkv2, pki2] = max(xcm2_pos);
        if pkv2 >= 4 * median(xcm2_pos)
            if pki2 > 1 && pki2 < numel(xcm2_pos)
                yA2 = xcm2_pos(pki2-1);  yB2 = xcm2_pos(pki2);  yC2 = xcm2_pos(pki2+1);
                d2 = 0.5 * (yA2 - yC2) / (yA2 - 2*yB2 + yC2);
                pk2 = lg2_pos(pki2) + d2;
            else
                pk2 = lg2_pos(pki2);
            end
            err_noDop(mc) = (pk2 - ref_delay_smp) * range_per_sample;
        end
    end

    % Compute RMSE for this SNR point
    valid_interp = err_interp(~isnan(err_interp));
    valid_raw    = err_raw   (~isnan(err_raw));
    valid_noDop  = err_noDop (~isnan(err_noDop));
    if ~isempty(valid_interp)
        range_RMSE_interp(sIdx) = sqrt(mean(valid_interp.^2));
        range_RMSE_raw   (sIdx) = sqrt(mean(valid_raw.^2));
    end
    if ~isempty(valid_noDop)
        range_RMSE_noDopComp(sIdx) = sqrt(mean(valid_noDop.^2));
    end
    detect_rate(sIdx) = detected / N_MC_sweep;
end

% ---- Plot ----
figure('Name','Ranging RMSE vs SNR','Position',[250 200 950 650],'Color','w');

yyaxis left        % Left Y axis: RMSE [m], log scale
semilogy(SNR_sweep_dB, range_RMSE_interp, 'b-', 'LineWidth', 2.5, ...
         'DisplayName', 'RMSE — Doppler comp + parabolic'); hold on;
semilogy(SNR_sweep_dB, range_RMSE_noDopComp, 'r-', 'LineWidth', 1.8, ...
         'DisplayName', sprintf('RMSE — NO Doppler comp (f_d=%.0f kHz)', doppler_test_Hz/1e3));
semilogy(SNR_sweep_dB, CRLB_dist, 'k--', 'LineWidth', 1.8, ...
         'DisplayName', sprintf('CRLB (B=%.1f MHz, L=%d)', BW_prs/1e6, L_prs));
yline(quant_floor_m, 'm:', 'Sample resolution', 'LineWidth', 1.5, 'HandleVisibility','off');
grid on; box on;
ylabel('Range RMSE (m)', 'FontWeight','bold','FontSize',11);
ylim([1e-2 1e5]);
ax = gca;  ax.YColor = 'k';

yyaxis right       % Right Y axis: detection rate [%]
plot(SNR_sweep_dB, detect_rate*100, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
     'DisplayName', 'Detection rate');
ylabel('Detection rate (%)', 'FontWeight','bold','FontSize',11);
ylim([0 110]);  ax.YColor = [0.4 0.4 0.4];

xlabel('SNR (dB)','FontWeight','bold','FontSize',11);
title('PRS Ranging Accuracy vs SNR (with Doppler)', 'FontSize',13,'FontWeight','bold');
legend('Location','southwest','FontSize',9);

fprintf('\n  BW_prs          : %.2f MHz   (full active pilot span)\n', BW_prs/1e6);
fprintf('  L_prs           : %d samples (matched-filter length)\n', L_prs);
fprintf('  Doppler test    : %.0f kHz\n',  doppler_test_Hz/1e3);
fprintf('  CRLB @ SNR=0 dB : %.4f m\n',   CRLB_dist(SNR_sweep_dB==0));
fprintf('  CRLB @ SNR=20dB : %.4f m\n',   CRLB_dist(SNR_sweep_dB==20));

% ----- export results into the shared context ------------------------
ctx.SNR_sweep_dB=SNR_sweep_dB; ctx.CRLB_dist=CRLB_dist;
ctx.range_RMSE_interp=range_RMSE_interp; ctx.range_RMSE_raw=range_RMSE_raw;
ctx.range_RMSE_noDopComp=range_RMSE_noDopComp; ctx.detect_rate=detect_rate;

end
