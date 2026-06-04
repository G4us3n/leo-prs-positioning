function ctx = stage10_montecarlo(cfg, ctx)
%STAGE10_MONTECARLO  Monte Carlo cross-correlation ranging (ALL visible satellites).
%
%   Pipeline stage. Reads CTX/CFG, runs MC ranging for every visible
%   satellite over all reference epochs, and returns updated CTX with 
%   detailed matrices (RMSE, SNR, Range, Elev, Doppler) for all links.

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

fprintf('\n>>> Section 10: MC Ranging with Doppler (ALL visible satellites) <<<\n');
N_MC       = cfg.mc.nTrials;        % Monte Carlo iterations
acc_RMSE_m = cfg.mc.accTarget_m;    % Accuracy target (meters)

for r = 1:nRef
    fprintf('\n=========================================================================================\n');
    fprintf(' DETAILED RANGING LOG | %s\n', ep(r).label);
    fprintf('-----------------------------------------------------------------------------------------\n');
    fprintf('%-10s | %-6s | %-9s | %-8s | %-10s | %-9s | %-8s | %-6s\n', ...
        'City','Sat','SNR(dB)','Elev','Range(km)','RMSE(m)','Doppler','STATUS');
    fprintf('-----------------------------------------------------------------------------------------\n');
    
    rIdx = refTimeIdx(r);
    
    % Initialize 2D matrices [Cities x Sats] to store ALL detailed link data
    ep(r).RMSE        = nan(numCities, nSats);
    ep(r).detectRate  = nan(numCities, nSats);
    ep(r).SNR_all     = nan(numCities, nSats);
    ep(r).Range_all   = nan(numCities, nSats);
    ep(r).Elev_all    = nan(numCities, nSats);
    ep(r).Doppler_all = nan(numCities, nSats);
    
    total_visible = 0; 
    nFail = 0;
    best_snr = -Inf;
    best_loc = '';
    best_rmse = NaN;
    
    for cIdx = 1:numCities
        cityName = cityNetwork{cIdx,1};
        
        for sIdx = 1:nSats
            elev = history_all_Elev(rIdx, cIdx, sIdx);
            
            % Skip satellites below the 5-degree horizon mask
            if isnan(elev) || elev < 5 
                continue; 
            end
            
            total_visible = total_visible + 1;
            
            % Extract physical parameters for this specific link
            SNRdB      = history_all_SNR(rIdx, cIdx, sIdx);
            true_d_m   = history_all_Range(rIdx, cIdx, sIdx) * 1000;
            doppler_Hz = history_all_Doppler(rIdx, cIdx, sIdx);
            if isnan(doppler_Hz), doppler_Hz = 0; end
            
            % Store the detailed parameters into the epoch struct
            ep(r).SNR_all(cIdx, sIdx)     = SNRdB;
            ep(r).Range_all(cIdx, sIdx)   = true_d_m;
            ep(r).Elev_all(cIdx, sIdx)    = elev;
            ep(r).Doppler_all(cIdx, sIdx) = doppler_Hz;
            
            % ---- Build the received signal ----
            ref_delay_smp = true_d_m / range_per_sample;    % True delay [samples]
            intRef  = floor(ref_delay_smp);                 % Integer part
            fracRef = ref_delay_smp - intRef;               % Fractional part
            
            L_mc = L_prs + intRef + 200;                    % Padded length
            Smc  = fft([prs_tx; zeros(L_mc - L_prs, 1)]);   % FFT of padded PRS
            kmc  = (0:L_mc-1).';
            kmc(kmc > L_mc/2) = kmc(kmc > L_mc/2) - L_mc;
            
            % Apply fractional delay in freq. domain
            sig_frac = ifft(Smc .* exp(-1j*2*pi*kmc*fracRef/L_mc));
            rx_t = circshift(sig_frac, intRef);   
            rx_t(1:intRef) = 0;                   
            
            % Apply Doppler frequency shift
            t_vec = (0:L_mc-1).' * Ts_prs;
            rx_t  = rx_t .* exp(1j*2*pi*doppler_Hz*t_vec);
            
            % Calculate noise power based on target SNR
            sig_pow = mean(abs(sig_frac).^2);
            n_pow   = sig_pow / 10^(SNRdB/10);
            
            err_smp  = nan(N_MC, 1);
            detected = 0;
            
            % ---- Monte Carlo Trials ----
            for mc = 1:N_MC
                noise_v = sqrt(n_pow/2) * (randn(L_mc,1) + 1j*randn(L_mc,1));
                rx_mc   = rx_t + noise_v;
                
                % Receiver applies Doppler compensation before correlation
                rx_mc_comp = rx_mc .* exp(-1j*2*pi*doppler_Hz*t_vec);
                [xcm, lgm] = xcorr(rx_mc_comp, prs_tx);
                xcm_pos = abs(xcm(lgm >= 0));
                lg_pos  = lgm(lgm >= 0);
                
                [pkv, pki] = max(xcm_pos);
                
                % Peak detection gate (must exceed threshold)
                if pkv >= cfg.mc.detectionGate * median(xcm_pos)     
                    % Sub-sample parabolic interpolation
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
            
            % Calculate RMSE
            valid_e = err_smp(~isnan(err_smp));
            if ~isempty(valid_e)
                curr_RMSE = sqrt(mean((valid_e * range_per_sample).^2));
            else
                curr_RMSE = NaN;
            end
            
            % Store results
            ep(r).RMSE(cIdx, sIdx)       = curr_RMSE;
            ep(r).detectRate(cIdx, sIdx) = detected / N_MC;
            
            if isnan(curr_RMSE) || curr_RMSE > acc_RMSE_m
                status = 'FAIL';  
                nFail = nFail + 1;
            else
                status = 'OK';
            end
            
            % Track the strongest overall link
            if SNRdB > best_snr
                best_snr  = SNRdB;
                best_loc  = sprintf('%s / %s', cityName, satNames{sIdx});
                best_rmse = curr_RMSE;
            end
            
            % Print detailed log ONLY for Milan (cIdx == 1) to prevent console spam
            if cIdx == 1
                fprintf('%-10s | %-6s | %+6.2f dB | %6.1f° | %7.2f km |  %7.2f m | %+.1fkHz | %-6s\n', ...
                    cityName, satNames{sIdx}, SNRdB, elev, true_d_m/1000, curr_RMSE, doppler_Hz/1e3, status);
            end
        end 
    end 
    
    fprintf('-----------------------------------------------------------------------------------------\n');
    fprintf(' Epoch Summary: %d visible links tested. Failing %dm threshold: %d links.\n', ...
            total_visible, acc_RMSE_m, nFail);
    fprintf(' Strongest overall link: %s (%.2f dB) | RMSE = %.2f m\n', best_loc, best_snr, best_rmse);

    % SNR-after-FFT comparison for Milan's visible links (paper §II-B, eq. 7)
    milan_snr_vec = ep(r).SNR_all(1, ~isnan(ep(r).SNR_all(1,:)));
    if ~isempty(milan_snr_vec)
        [fft_snr_vec, ~, ~, ~] = snr_after_fft(milan_snr_vec, L_prs, BW_prs, c);
        fprintf(' SNR-after-FFT Milan (+%.0f dB gain): per-sample [%.1f, %.1f] dB  =>  post-FFT [%.1f, %.1f] dB\n', ...
            10*log10(L_prs), min(milan_snr_vec), max(milan_snr_vec), ...
            min(fft_snr_vec), max(fft_snr_vec));
    end
end

% Export updated struct array back into context
ctx.ep = ep;  
end