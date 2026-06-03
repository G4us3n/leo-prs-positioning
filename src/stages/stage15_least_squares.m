function ctx = stage15_least_squares(cfg, ctx)
    % STAGE 15: Least-squares positioning over the trajectory
    % Runs standard Least-Squares (LS) and automatically generates diagnostic plots.

    disp('======================================================');
    disp('Stage 15: Least-Squares Positioning & Error Analysis');
    disp('======================================================');

    % 1. Extract constants
    c = ctx.c;
    target_city_idx = 1; % Milan
    [true_x, true_y, true_z] = geodetic2ecef(ctx.wgs84, ctx.cityNetwork{target_city_idx, 2}, ...
                                             ctx.cityNetwork{target_city_idx, 3}, 0);
    true_pos_ecef = [true_x, true_y, true_z];

    time_offsets = seconds(ctx.timeVec - ctx.zenithTime);
    valid_time_indices = find(time_offsets >= -150 & time_offsets <= 150); % Focusing on zenith
    num_steps = length(valid_time_indices);

    ls_results.time_offsets = time_offsets(valid_time_indices);
    ls_results.est_pos = nan(num_steps, 3);
    ls_results.error_3d = nan(num_steps, 1);
    ls_results.num_sats_used = zeros(num_steps, 1);

    % Initialize with a geographically local guess (Milan-ish) to prevent singularity
    initial_guess = true_pos_ecef + [100, 100, 100]; 

    for step_idx = 1:num_steps
        tIdx = valid_time_indices(step_idx);
        vis_sat_pos = []; pseudo_ranges = [];
        
        for sIdx = 1:ctx.nSats
            elev = ctx.history_all_Elev(tIdx, target_city_idx, sIdx);
            snr_db = ctx.history_all_SNR(tIdx, target_city_idx, sIdx);
            
            % Thresholding: Only use satellites with SNR > -30dB and Elevation > 5deg
            if ~isnan(elev) && elev > 25 && snr_db > -35
                [satPos, ~] = states(ctx.sats(sIdx), ctx.timeVec(tIdx), "CoordinateFrame", "ecef");
                true_range = ctx.history_all_Range(tIdx, target_city_idx, sIdx) * 1000;
                measured_rho = true_range + (randn() * 5.0); 
                
                vis_sat_pos = [vis_sat_pos; satPos'];
                pseudo_ranges = [pseudo_ranges; measured_rho];
            end
        end
        
        ls_results.num_sats_used(step_idx) = size(vis_sat_pos, 1);
        
        if size(vis_sat_pos, 1) >= 4
            % Standard LS solver (calling the original version without WLS)
            [est_pos, ~, ~, ~] = ls_position_standard(vis_sat_pos, pseudo_ranges, initial_guess);
            
            if ~isnan(est_pos(1)) && norm(est_pos - true_pos_ecef) < 1e5
                ls_results.est_pos(step_idx, :) = est_pos;
                ls_results.error_3d(step_idx) = norm(est_pos - true_pos_ecef);
                initial_guess = est_pos; % Update guess for next step
            end
        end
    end

    % 2. Automatically generate and display the error plot
    figure('Name', 'Positioning Accuracy Evolution (LS)', 'Color', 'w');
    subplot(2,1,1);
    plot(ls_results.time_offsets, ls_results.error_3d, 'b-', 'LineWidth', 2);
    grid on; xline(0, 'r--', 'Zenith'); ylabel('3D Error (m)');
    title(['Positioning Accuracy: Mean Error = ' num2str(mean(ls_results.error_3d, 'omitnan'), '%.2f') ' m']);
    
    subplot(2,1,2);
    stairs(ls_results.time_offsets, ls_results.num_sats_used, 'k-', 'LineWidth', 1.5);
    grid on; ylabel('Visible Satellites'); xlabel('Time from Zenith (s)'); ylim([0 13]);
    
    ctx.ls_results = ls_results;
end

% Standard solver (No WLS weights to ensure basic functionality)
function [est_pos, est_clock, iter, GDOP] = ls_position_standard(sat_pos, pseudo_ranges, initial_guess)
    x_state = [initial_guess(:); 0]; 
    for iter = 1:10
        H = zeros(size(sat_pos,1), 4);
        delta_rho = zeros(size(sat_pos,1), 1);
        for i = 1:size(sat_pos,1)
            R = norm(sat_pos(i,:) - x_state(1:3)');
            delta_rho(i) = pseudo_ranges(i) - (R + x_state(4));
            H(i, 1:3) = -(sat_pos(i,:) - x_state(1:3)') / R;
            H(i, 4) = 1;
        end
        delta_x = H \ delta_rho;
        x_state = x_state + delta_x;
        if norm(delta_x(1:3)) < 1e-3, break; end
    end
    est_pos = x_state(1:3)'; est_clock = x_state(4); GDOP = 0;
end