function ctx = stage13_dop(cfg, ctx)
%STAGE13_DOP  Geometric Dilution of Precision (HDOP/VDOP/PDOP) at each epoch.
%
%   Pipeline stage. Reads CTX/CFG, runs one section, returns updated CTX.

sats=ctx.sats; nSats=ctx.nSats; nPlanes=ctx.nPlanes; numCities=ctx.numCities;
cityNetwork=ctx.cityNetwork; wgs84=ctx.wgs84; timeVec=ctx.timeVec;
history_all_Elev=ctx.history_all_Elev; ep=ctx.ep; nRef=ctx.nRef;
refTimeIdx=ctx.refTimeIdx; allTags=ctx.allTags;

% ----- original section body (unchanged physics) ---------------------
fprintf('\n>>> Section 11: Geometric DOP Analysis <<<\n');

isCoplanar = (nPlanes == 1);
if isCoplanar
    fprintf('  Single plane → computing 2D HDOP only (E,N,clock)\n\n');
else
    fprintf('  %d orbital planes → computing full 3D DOP (E,N,U,clock)\n\n', nPlanes);
end

fprintf('------------------------------------------------------------------------\n');
if isCoplanar
    fprintf('%-10s | %-16s | %-8s | %-8s | %-6s\n', 'City','Epoch','HDOP','RCOND','nSats');
else
    fprintf('%-10s | %-16s | %-8s | %-8s | %-8s | %-6s\n', 'City','Epoch','HDOP','VDOP','PDOP','nSats');
end
fprintf('------------------------------------------------------------------------\n');

for r = 1:nRef
    rIdx = refTimeIdx(r);

    for cIdx = 1:numCities
        if ~ep(r).visMask(cIdx), continue; end
        cLat = cityNetwork{cIdx,2};  cLon = cityNetwork{cIdx,3};
        [cX, cY, cZ] = geodetic2ecef(wgs84, cLat, cLon, 0);

        % Build ENU (East-North-Up) frame at the ground station
        az_sats = [];  el_sats = [];
        R_gs  = norm([cX, cY, cZ]);
        up    = [cX; cY; cZ] / R_gs;
        east  = [-cY; cX; 0];  east = east / norm(east);
        north = cross(up, east);

        % Compute azimuth and elevation to each visible satellite
        for sIdx = 1:nSats
            if isnan(history_all_Elev(rIdx, cIdx, sIdx)), continue; end
            [satPos, ~] = states(sats(sIdx), timeVec(rIdx), "CoordinateFrame","ecef");
            V      = satPos - [cX; cY; cZ];
            V_unit = V / norm(V);
            az_sats(end+1) = atan2d(dot(V_unit, east), dot(V_unit, north));
            el_sats(end+1) = asind(dot(V_unit, up));
        end

        nVisSats = numel(az_sats);
        if nVisSats < 2, continue; end

        if isCoplanar
            % 2D: 3 unknowns (E, N, Clock)
            H2 = zeros(nVisSats, 3);
            for s = 1:nVisSats
                H2(s,1) = cosd(el_sats(s)) * sind(az_sats(s));
                H2(s,2) = cosd(el_sats(s)) * cosd(az_sats(s));
                H2(s,3) = 1;
            end
            G  = H2' * H2;
            rc = rcond(G);
            if rc > 1e-12
                Q2   = inv(G);
                HDOP = sqrt(Q2(1,1) + Q2(2,2));
                fprintf('%-10s | %-16s | %6.2f   | %.2e | %d\n', ...
                    cityNetwork{cIdx,1}, allTags{r}, HDOP, rc, nVisSats);
            else
                fprintf('%-10s | %-16s | ill-cond | %.2e | %d\n', ...
                    cityNetwork{cIdx,1}, allTags{r}, rc, nVisSats);
            end
        else
            % 3D: 4 unknowns (E, N, U, Clock) — need ≥ 4 visible sats
            H4 = zeros(nVisSats, 4);
            for s = 1:nVisSats
                H4(s,1) = cosd(el_sats(s)) * sind(az_sats(s));   % East
                H4(s,2) = cosd(el_sats(s)) * cosd(az_sats(s));   % North
                H4(s,3) = sind(el_sats(s));                        % Up
                H4(s,4) = 1;                                        % Clock bias
            end
            G  = H4' * H4;
            rc = rcond(G);
            if rc > 1e-12 && nVisSats >= 4
                Q4   = inv(G);
                HDOP = sqrt(Q4(1,1) + Q4(2,2));
                VDOP = sqrt(Q4(3,3));
                PDOP = sqrt(Q4(1,1) + Q4(2,2) + Q4(3,3));
                fprintf('%-10s | %-16s | %6.2f   | %6.2f   | %6.2f   | %d\n', ...
                    cityNetwork{cIdx,1}, allTags{r}, HDOP, VDOP, PDOP, nVisSats);
            else
                fprintf('%-10s | %-16s | need 4+ sats on diverse planes  | %d\n', ...
                    cityNetwork{cIdx,1}, allTags{r}, nVisSats);
            end
        end
    end
end
fprintf('------------------------------------------------------------------------\n');
if isCoplanar
    fprintf('Note: Single-plane constellation → only 2D horizontal DOP is computable.\n');
    fprintf('      Add more orbital planes to enable full 3D DOP.\n');
end

% ----- export results into the shared context ------------------------
% (DOP stage prints a table; no exports)

end
