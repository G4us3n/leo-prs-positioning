function ctx = stage05_flyover(cfg, ctx)
%STAGE05_FLYOVER  Time-series flyover: per-link elevation, range, power, SNR, Doppler.
%
%   Pipeline stage. Reads CTX/CFG, runs one section, returns updated CTX.

sats=ctx.sats; gs=ctx.gs; nSats=ctx.nSats; nPlanes=ctx.nPlanes;
nSatsPerPlane=ctx.nSatsPerPlane; numCities=ctx.numCities; cityNetwork=ctx.cityNetwork;
wgs84=ctx.wgs84; zenithTime=ctx.zenithTime; satNames=ctx.satNames;
satPlaneIdx=ctx.satPlaneIdx; freq_Spare=ctx.freq_Spare; c=ctx.c; lambda=ctx.lambda;
pat3D=ctx.pat3D; AzMat=ctx.AzMat; ElMat=ctx.ElMat;
noise_floor_dBm=ctx.noise_floor_dBm; noise_density_dBm=ctx.noise_density_dBm;
T_ref=ctx.T_ref; cfg_atm=ctx.cfg_atm;

% ----- original section body (unchanged physics) ---------------------
disp('======================================================');
fprintf('Initiating %d-Sat Time-Series Flyover Analysis (-3 min to +3 min)...\n', nSats);
disp('======================================================');

startTime  = zenithTime - minutes(3);
endTime    = zenithTime + minutes(3);
timeVec    = startTime : cfg.flyoverStep : endTime;   % Every 2 seconds
numEpochs  = length(timeVec);                     % 181 epochs
relTimeMin = seconds(timeVec - zenithTime) / 60;  % Relative time in minutes

% Pre-allocate 3-D arrays:  (epochs × cities × satellites)
history_all_Elev    = nan(numEpochs, numCities, nSats);   % Elevation [°]
history_all_Range   = nan(numEpochs, numCities, nSats);   % Slant range [km]
history_all_Power   = nan(numEpochs, numCities, nSats);   % Rx power [dBm]
history_all_SNR     = nan(numEpochs, numCities, nSats);   % SNR [dB]
history_all_Doppler = nan(numEpochs, numCities, nSats);   % Doppler shift [Hz]

% Satellite ground tracks (sub-satellite lat/lon over time)
sat_Trace_Lat = nan(numEpochs, nSats);
sat_Trace_Lon = nan(numEpochs, nSats);

for tIdx = 1:numEpochs
    currentTime = timeVec(tIdx);

    for sIdx = 1:nSats
        [satPos, satVel] = states(sats(sIdx), currentTime, "CoordinateFrame","ecef");
        [satLat, satLon, satAlt] = ecef2geodetic(wgs84, satPos(1), satPos(2), satPos(3));
        sat_Trace_Lat(tIdx, sIdx) = satLat;
        sat_Trace_Lon(tIdx, sIdx) = satLon;

        % Antenna frame (identical to Section 4)
        [nadirX, nadirY, nadirZ] = geodetic2ecef(wgs84, satLat, satLon, 0);
        X_ant = [nadirX; nadirY; nadirZ] - satPos;  X_ant = X_ant / norm(X_ant);
        Y_ant = cross(X_ant, satVel);                Y_ant = Y_ant / norm(Y_ant);
        Z_ant = cross(X_ant, Y_ant);

        for cIdx = 1:numCities
            cLat = cityNetwork{cIdx,2};  cLon = cityNetwork{cIdx,3};
            [cX, cY, cZ] = geodetic2ecef(wgs84, cLat, cLon, 0);

            % Vector from satellite to ground station
            V_x = cX - satPos(1);  V_y = cY - satPos(2);  V_z = cZ - satPos(3);
            SlantRange = norm([V_x, V_y, V_z]);
            V_ux = V_x/SlantRange;  V_uy = V_y/SlantRange;  V_uz = V_z/SlantRange;

            % Antenna-frame angles → gain look-up
            x_proj   = X_ant(1)*V_ux + X_ant(2)*V_uy + X_ant(3)*V_uz;
            y_proj   = Y_ant(1)*V_ux + Y_ant(2)*V_uy + Y_ant(3)*V_uz;
            z_proj   = Z_ant(1)*V_ux + Z_ant(2)*V_uy + Z_ant(3)*V_uz;
            az_local = atan2d(y_proj, x_proj);
            el_local = asind(z_proj);

            % Ground-station elevation angle
            R_norm    = norm([cX, cY, cZ]);
            Z_x = cX/R_norm;  Z_y = cY/R_norm;  Z_z = cZ/R_norm;
            sin_el_gs = (-V_x*Z_x - V_y*Z_y - V_z*Z_z) / SlantRange;
            el_gs_deg = asind(sin_el_gs);

            % Doppler shift
            %   LOS = line-of-sight unit vector from GS toward satellite
            %   v_radial = satellite velocity projected onto LOS
            %   Positive v_radial = satellite approaching → positive Doppler
            los_unit   = -[V_x, V_y, V_z] / SlantRange;
            v_radial   = dot(satVel, los_unit);
            doppler_Hz = (v_radial / c) * freq_Spare;

            if el_gs_deg >= 0    % Satellite above horizon
                % Atmospheric losses (ITU-R P.618 with gaspl fallback)
                el_clipped = max(el_gs_deg, 5);
                try
                    [atmLoss_total, ~, ~, ~, ~] = p618PropagationLosses( ...
                        freq_Spare, el_clipped, cLat, cLon, ...
                        cfg_atm.TotalColumnWaterVapourDensity, ...
                        cfg_atm.AnnualExceedance, ...
                        'Polarization', cfg_atm.Polarization);
                    atmLoss = atmLoss_total;
                catch
                    T_c = T_ref - 273.15;
                    atmLoss = gaspl(SlantRange, freq_Spare, T_c, 101.325, 7.5);
                end

                gain         = interp2(AzMat, ElMat, pat3D, az_local, el_local, 'linear', -60);
                prsPower_dBm = 10*log10(cfg.txPower_W) + 30;
                FSPL         = fspl(SlantRange, lambda);
                rxGain_dBi   = 10*log10(cfg.rxDishEff * (pi * cfg.rxDishDiam_m / lambda)^2);
                rxPwr        = prsPower_dBm + gain + rxGain_dBi - FSPL - atmLoss;

                history_all_Elev(tIdx, cIdx, sIdx)    = el_gs_deg;
                history_all_Range(tIdx, cIdx, sIdx)   = SlantRange / 1000;  % → km
                history_all_Power(tIdx, cIdx, sIdx)   = rxPwr;
                history_all_SNR(tIdx, cIdx, sIdx)     = rxPwr - noise_floor_dBm;
                history_all_Doppler(tIdx, cIdx, sIdx) = doppler_Hz;
            end
            % Satellite below horizon → entries stay NaN (= invisible)
        end
    end
end

% ---- Macro-diversity: best satellite per city per epoch ----
history_Elev  = nan(numEpochs, numCities);
history_Range = nan(numEpochs, numCities);
history_Power = nan(numEpochs, numCities);
history_SNR   = nan(numEpochs, numCities);

for tIdx = 1:numEpochs
    for cIdx = 1:numCities
        [max_snr, best_sIdx] = max(squeeze(history_all_SNR(tIdx, cIdx, :)));
        if ~isnan(max_snr)
            history_SNR(tIdx, cIdx)   = max_snr;
            history_Power(tIdx, cIdx) = history_all_Power(tIdx, cIdx, best_sIdx);
            history_Range(tIdx, cIdx) = history_all_Range(tIdx, cIdx, best_sIdx);
            history_Elev(tIdx, cIdx)  = history_all_Elev(tIdx, cIdx, best_sIdx);
        end
    end
end

% ----- export results into the shared context ------------------------
ctx.timeVec=timeVec; ctx.numEpochs=numEpochs; ctx.relTimeMin=relTimeMin;
ctx.history_all_Elev=history_all_Elev; ctx.history_all_Range=history_all_Range;
ctx.history_all_Power=history_all_Power; ctx.history_all_SNR=history_all_SNR;
ctx.history_all_Doppler=history_all_Doppler;
ctx.sat_Trace_Lat=sat_Trace_Lat; ctx.sat_Trace_Lon=sat_Trace_Lon;
ctx.history_Elev=history_Elev; ctx.history_Range=history_Range;
ctx.history_Power=history_Power; ctx.history_SNR=history_SNR;

end
