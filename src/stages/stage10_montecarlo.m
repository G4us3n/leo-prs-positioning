function ctx = stage10_montecarlo(cfg, ctx)
%STAGE10_MONTECARLO  Monte Carlo cross-correlation ranging (best satellite per city).
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
fprintf('\n>>> Section 8: MC Ranging with Doppler (best-sat per city) <<<\n');

N_MC       = cfg.mc.nTrials;   % Monte Carlo iterations
acc_RMSE_m = cfg.mc.accTarget_m;    % Accuracy target: 30 metres

for r = 1:nRef
    fprintf('\n=========================================================================================\n');
    fprintf(' RANGING QUALITY | %s\n', ep(r).label);
    fprintf('-----------------------------------------------------------------------------------------\n');
    fprintf('%-10s | %-7s | %-9s | %-8s | %-10s | %-9s | %-8s | %-6s\n', ...
        'City','BestSat','SNR(dB)','Elev','Range(km)','RMSE(m)','Doppler','STATUS');
    fprintf('-----------------------------------------------------------------------------------------\n');

    nFail = 0;  best_snr = -Inf;  best_loc = '';  best_rmse = NaN;
    rIdx = refTimeIdx(r);

    for cIdx = 1:numCities
        if ~ep(r).visMask(cIdx), continue; end

        cityName = cityNetwork{cIdx,1};
        bSat     = ep(r).bestSat(cIdx);
        SNRdB    = ep(r).SNR(cIdx);
        elev     = ep(r).Elev(cIdx);
        true_d_m = ep(r).Range(cIdx);

        doppler_Hz = history_all_Doppler(rIdx, cIdx, bSat);
        if isnan(doppler_Hz), doppler_Hz = 0; end

        % ---- Build the received signal ----
        ref_delay_smp = true_d_m / range_per_sample;   % True delay [samples]
        intRef  = floor(ref_delay_smp);                 % Integer part
        fracRef = ref_delay_smp - intRef;               % Fractional part

        L_mc = L_prs + intRef + 200;                    % Padded length
        Smc  = fft([prs_tx; zeros(L_mc - L_prs, 1)]);  % FFT of padded PRS
        kmc  = (0:L_mc-1).';
        kmc(kmc > L_mc/2) = kmc(kmc > L_mc/2) - L_mc;

        % Fractional delay in freq. domain
        sig_frac = ifft(Smc .* exp(-1j*2*pi*kmc*fracRef/L_mc));
        rx_t = circshift(sig_frac, intRef);   % Integer delay
        rx_t(1:intRef) = 0;                   % Zero the "wrapped" portion

        % Apply Doppler:  rx(t) × e^{j 2π f_d t}
        t_vec = (0:L_mc-1).' * Ts_prs;
        rx_t  = rx_t .* exp(1j*2*pi*doppler_Hz*t_vec);

        % Noise level for target SNR (per-sample convention)
        sig_pow = mean(abs(sig_frac).^2);
        n_pow   = sig_pow / 10^(SNRdB/10);

        err_smp  = nan(N_MC, 1);
        detected = 0;

        for mc = 1:N_MC
            noise_v = sqrt(n_pow/2) * (randn(L_mc,1) + 1j*randn(L_mc,1));
            rx_mc   = rx_t + noise_v;

            % Doppler compensation
            rx_mc_comp = rx_mc .* exp(-1j*2*pi*doppler_Hz*t_vec);

            [xcm, lgm] = xcorr(rx_mc_comp, prs_tx);
            xcm_pos = abs(xcm(lgm >= 0));
            lg_pos  = lgm(lgm >= 0);
            [pkv, pki] = max(xcm_pos);

            if pkv >= cfg.mc.detectionGate*median(xcm_pos)     % Detection gate
                % Parabolic interpolation
                if pki > 1 && pki < numel(xcm_pos)
                    yA = xcm_pos(pki-1);  yB = xcm_pos(pki);  yC = xcm_pos(pki+1);
                    d  = 0.5*(yA - yC) / (yA - 2*yB + yC);
                    pk_l = lg_pos(pki) + d;
                else
                    pk_l = lg_pos(pki);
                end
                err_smp(mc) = pk_l - ref_delay_smp;
                detected    = detected + 1;
            end
        end

        valid_e = err_smp(~isnan(err_smp));
        if ~isempty(valid_e)
            curr_RMSE = sqrt(mean((valid_e * range_per_sample).^2));
        else
            curr_RMSE = NaN;
        end

        ep(r).RMSE(cIdx)       = curr_RMSE;
        ep(r).detectRate(cIdx) = detected / N_MC;

        if isnan(curr_RMSE) || curr_RMSE > acc_RMSE_m
            status = 'FAIL';  nFail = nFail + 1;
        else
            status = 'OK';
        end
        if SNRdB > best_snr
            best_snr  = SNRdB;
            best_loc  = sprintf('%s / %s', cityName, satNames{bSat});
            best_rmse = curr_RMSE;
        end

        fprintf('%-10s | %-6s | %+6.2f dB | %6.1f° | %7.2f km |  %7.2f m | %+.1fkHz | %-6s\n', ...
            cityName, satNames{bSat}, SNRdB, elev, true_d_m/1000, curr_RMSE, doppler_Hz/1e3, status);
    end

    fprintf('-----------------------------------------------------------------------------------------\n');
    fprintf(' Strongest: %s (%.2f dB) | RMSE=%.2f m | Failing %dm: %d/%d\n', ...
        best_loc, best_snr, best_rmse, acc_RMSE_m, nFail, sum(ep(r).visMask));
end

% ----- export results into the shared context ------------------------
ctx.ep=ep;  % RMSE and detectRate fields now populated

end
