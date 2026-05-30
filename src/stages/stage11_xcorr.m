function ctx = stage11_xcorr(cfg, ctx)
%STAGE11_XCORR  Cross-correlation visualisation for representative links.
%
%   Pipeline stage. Reads CTX/CFG, runs one section, returns updated CTX.

sats=ctx.sats; nSats=ctx.nSats; numCities=ctx.numCities; cityNetwork=ctx.cityNetwork;
wgs84=ctx.wgs84; zenithTime=ctx.zenithTime; satNames=ctx.satNames;
freq_Spare=ctx.freq_Spare; c=ctx.c; lambda=ctx.lambda;
pat3D=ctx.pat3D; AzMat=ctx.AzMat; ElMat=ctx.ElMat;
noise_floor_dBm=ctx.noise_floor_dBm; T_ref=ctx.T_ref; cfg_atm=ctx.cfg_atm;
timeVec=ctx.timeVec; relTimeMin=ctx.relTimeMin;
history_all_Elev=ctx.history_all_Elev; history_all_Range=ctx.history_all_Range;
history_all_SNR=ctx.history_all_SNR; history_all_Doppler=ctx.history_all_Doppler;
prs_tx=ctx.prs_tx; Ts_prs=ctx.Ts_prs; range_per_sample=ctx.range_per_sample;
L_prs=ctx.L_prs; BW_prs=ctx.BW_prs;
ep=ctx.ep; nRef=ctx.nRef; refTimeIdx=ctx.refTimeIdx; allTags=ctx.allTags;

% ----- original section body (unchanged physics) ---------------------
for r = 1:nRef
    visIdx = find(ep(r).visMask);
    if numel(visIdx) < 3, continue; end

    [~, ord] = sort(ep(r).SNR(visIdx), 'descend');
    pickIdx  = visIdx(ord(round(linspace(1, numel(ord), 3))));
    pickLbl  = {'Strongest', 'Median', 'Weakest'};

    figure('Name', ['PRS ToA - ' ep(r).label], ...
           'Position',[60+30*r, 60+30*r, 1300, 720], 'Color', 'w');
    sgtitle(sprintf('PRS Cross-Correlation @ %s', ep(r).label), 'FontWeight','bold');

    for k = 1:3
        cIdx  = pickIdx(k);
        bSat  = ep(r).bestSat(cIdx);
        SNRdB = ep(r).SNR(cIdx);
        trueR = ep(r).Range(cIdx);
        elv   = ep(r).Elev(cIdx);

        trueDsmp = trueR / range_per_sample;
        intD = floor(trueDsmp);  fracD = trueDsmp - intD;

        L_total = L_prs + intD + 100;
        Sf = fft([prs_tx; zeros(L_total - L_prs, 1)]);
        kv = (0:L_total-1).';  kv(kv > L_total/2) = kv(kv > L_total/2) - L_total;
        sig_fd = ifft(Sf .* exp(-1j*2*pi*kv*fracD/L_total));
        rx_t = circshift(sig_fd, intD);  rx_t(1:intD) = 0;

        nPow = mean(abs(sig_fd).^2) / 10^(SNRdB/10);
        noise_v = sqrt(nPow/2) * (randn(size(rx_t)) + 1j*randn(size(rx_t)));
        rx_sig  = rx_t + noise_v;

        [xc, lags] = xcorr(rx_sig, prs_tx);
        xc_pos = abs(xc(lags>=0));  lg_pos = lags(lags>=0);
        [~, pi_] = max(xc_pos);
        if pi_>1 && pi_<numel(xc_pos)
            yA = xc_pos(pi_-1);  yB = xc_pos(pi_);  yC = xc_pos(pi_+1);
            d  = 0.5*(yA - yC) / (yA - 2*yB + yC);
            pk_l = lg_pos(pi_) + d;
        else
            pk_l = lg_pos(pi_);
        end
        est_R = pk_l * range_per_sample;
        err_R = est_R - trueR;

        delay_us = lg_pos * Ts_prs * 1e6;
        range_km = lg_pos * Ts_prs * c / 1e3;

        subplot(2,3,k);
        plot(delay_us, 20*log10(xc_pos/max(xc_pos)+eps), 'LineWidth', 1); hold on;
        xline(trueDsmp*Ts_prs*1e6, 'g--', 'True', 'LineWidth', 1.5);
        xline(pk_l*Ts_prs*1e6,     'r:',  'Est',  'LineWidth', 1.5);
        grid on; box on; ylim([-40 5]);
        xlim(trueDsmp*Ts_prs*1e6 + [-5 5]);
        xlabel('Delay (\mus)'); ylabel('|XCorr| (dB)');
        title(sprintf('%s | %s | %s\nSNR=%.1f dB | Elev=%.1f°', ...
              cityNetwork{cIdx,1}, satNames{bSat}, pickLbl{k}, SNRdB, elv), 'FontWeight','bold');

        subplot(2,3,k+3);
        plot(range_km - trueR/1000, xc_pos/max(xc_pos), 'LineWidth', 1); hold on;
        xline(0,                         'g--', 'True', 'LineWidth', 1.5);
        xline(est_R/1000 - trueR/1000,  'r:',  'Est',  'LineWidth', 1.5);
        grid on; box on; xlim([-0.3 0.3]);
        xlabel('Range error (km)'); ylabel('Norm correlation');
        title(sprintf('Range error = %.2f m', err_R));
    end
end

% ----- export results into the shared context ------------------------
% (cross-correlation stage produces figures only)

end
