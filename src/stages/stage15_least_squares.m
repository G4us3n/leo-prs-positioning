function ctx = stage15_least_squares(cfg, ctx)
% STAGE15_LEAST_SQUARES  LS positioning over the trajectory for every city.
%
%   For each city in the network, runs a standard iterative Least-Squares
%   (LS) solver at each epoch in the ±150 s window around zenith, then
%   calls plot_ls_accuracy_evolution to display 3D error, horizontal error,
%   and satellite count per city.

disp('======================================================');
disp('Stage 15: Least-Squares Positioning & Error Analysis');
disp('======================================================');

time_offsets       = seconds(ctx.timeVec - ctx.zenithTime);
valid_time_indices = find(time_offsets >= -150 & time_offsets <= 150);
num_steps          = length(valid_time_indices);

all_ls_results = cell(ctx.numCities, 1);

for city_idx = 1:ctx.numCities
    city_name = ctx.cityNetwork{city_idx, 1};
    fprintf('  LS positioning: %s ...\n', city_name);

    [true_x, true_y, true_z] = geodetic2ecef(ctx.wgs84, ...
        ctx.cityNetwork{city_idx, 2}, ctx.cityNetwork{city_idx, 3}, 0);
    true_pos_ecef = [true_x, true_y, true_z];
    initial_guess = true_pos_ecef + [100, 100, 100];

    r.time_offsets     = time_offsets(valid_time_indices);
    r.est_pos          = nan(num_steps, 3);
    r.error_horizontal = nan(num_steps, 1);
    r.error_3d         = nan(num_steps, 1);
    r.num_sats_used    = zeros(num_steps, 1);

    for step_idx = 1:num_steps
        tIdx = valid_time_indices(step_idx);
        vis_sat_pos = []; pseudo_ranges = [];

        for sIdx = 1:ctx.nSats
            elev   = ctx.history_all_Elev(tIdx, city_idx, sIdx);
            snr_db = ctx.history_all_SNR(tIdx, city_idx, sIdx);

            if ~isnan(elev) && elev > 25 && snr_db > -35
                [satPos, ~] = states(ctx.sats(sIdx), ctx.timeVec(tIdx), ...
                                     "CoordinateFrame", "ecef");
                true_range   = ctx.history_all_Range(tIdx, city_idx, sIdx) * 1000;
                measured_rho = true_range + randn() * 5.0;

                vis_sat_pos   = [vis_sat_pos;   satPos'];
                pseudo_ranges = [pseudo_ranges; measured_rho];
            end
        end

        r.num_sats_used(step_idx) = size(vis_sat_pos, 1);

        if size(vis_sat_pos, 1) >= 4
            [est_pos, ~, ~, ~] = ls_position_standard(vis_sat_pos, pseudo_ranges, initial_guess);

            if ~isnan(est_pos(1)) && norm(est_pos - true_pos_ecef) < 1e5
                r.est_pos(step_idx, :)      = est_pos;
                pos_error                   = est_pos - true_pos_ecef;
                r.error_3d(step_idx)        = norm(pos_error);
                r.error_horizontal(step_idx) = norm(pos_error(1:2));
                initial_guess               = est_pos;
            end
        end
    end

    all_ls_results{city_idx} = r;
    plot_ls_accuracy_evolution(r, city_name);
end

ctx.ls_results     = all_ls_results{1};   % Milan kept as primary result (backward compat)
ctx.all_ls_results = all_ls_results;
end

% -------------------------------------------------------------------------
function [est_pos, est_clock, iter, GDOP] = ls_position_standard(sat_pos, pseudo_ranges, initial_guess)
% Standard iterative LS solver (no WLS weights).
x_state = [initial_guess(:); 0];
for iter = 1:10
    H = zeros(size(sat_pos, 1), 4);
    delta_rho = zeros(size(sat_pos, 1), 1);
    for i = 1:size(sat_pos, 1)
        R = norm(sat_pos(i,:) - x_state(1:3)');
        delta_rho(i) = pseudo_ranges(i) - (R + x_state(4));
        H(i, 1:3)    = -(sat_pos(i,:) - x_state(1:3)') / R;
        H(i, 4)      = 1;
    end
    delta_x = H \ delta_rho;
    x_state = x_state + delta_x;
    if norm(delta_x(1:3)) < 1e-3, break; end
end
est_pos   = x_state(1:3)';
est_clock = x_state(4);
GDOP      = 0;
end

% -------------------------------------------------------------------------
function plot_ls_accuracy_evolution(ls_results, city_name)
% PLOT_LS_ACCURACY_EVOLUTION  3-panel accuracy plot for one city.
t     = ls_results.time_offsets;
err3d = ls_results.error_3d;
errH  = ls_results.error_horizontal;
nSat  = ls_results.num_sats_used;

rmse3d = sqrt(mean(err3d.^2, 'omitnan'));
mean3d = mean(err3d,         'omitnan');
r95_3d = prctile(err3d, 95);

figure('Name', sprintf('LS Accuracy — %s', city_name), 'Color', 'w');

subplot(3,1,1);
plot(t, err3d, 'b-', 'LineWidth', 2);
grid on; xline(0, 'r--', 'Zenith');
ylabel('3D Error (m)');
title(sprintf('%s  |  Mean = %.1f m   RMSE = %.1f m   R95 = %.1f m', ...
    city_name, mean3d, rmse3d, r95_3d));

subplot(3,1,2);
plot(t, errH, 'm-', 'LineWidth', 2);
grid on; xline(0, 'r--', 'Zenith');
ylabel('Horizontal Error (m)');
title('Horizontal Positioning Error');

subplot(3,1,3);
stairs(t, nSat, 'k-', 'LineWidth', 1.5);
grid on; xline(0, 'r--', 'Zenith');
ylabel('Satellites Used'); xlabel('Time from Zenith (s)');
ylim([0 max(nSat) + 1]);
end
