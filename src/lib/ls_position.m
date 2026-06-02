function [est_pos, est_clock, iter, GDOP] = ls_position(sat_pos, pseudo_ranges, initial_guess, snr_db, prs_bw)
    % LS_POSITION: Weighted Least Squares (WLS) solver for 3D positioning
    % Inputs:
    %   sat_pos       - N x 3 matrix of visible satellite ECEF coordinates
    %   pseudo_ranges - N x 1 vector of measured pseudoranges (meters)
    %   initial_guess - 1 x 3 vector of starting position guess
    %   snr_db        - N x 1 vector of SNR for each link (dB)
    %   prs_bw        - PRS active bandwidth (Hz)
    % Outputs:
    %   est_pos, est_clock, iter, GDOP
    
    max_iter = 15;
    tol = 1e-3; 
    
    % Constants for CRLB noise model
    c = 299792458; 
    
    x_state = [initial_guess(:); 0]; 
    num_sats = length(pseudo_ranges);
    
    if num_sats < 4, [est_pos, est_clock, iter, GDOP] = deal([NaN,NaN,NaN], NaN, 0, NaN); return; end
    
    for iter = 1:max_iter
        H = zeros(num_sats, 4);
        delta_rho = zeros(num_sats, 1);
        W = zeros(num_sats, num_sats); % Weighting matrix
        
        rec_pos = x_state(1:3)'; 
        clock_bias = x_state(4);
        
        for i = 1:num_sats
            dx = sat_pos(i, 1) - rec_pos(1);
            dy = sat_pos(i, 2) - rec_pos(2);
            dz = sat_pos(i, 3) - rec_pos(3);
            R_guess = sqrt(dx^2 + dy^2 + dz^2);
            
            delta_rho(i) = pseudo_ranges(i) - (R_guess + clock_bias);
            
            H(i, 1:3) = -[dx, dy, dz] / R_guess;
            H(i, 4) = 1;
            
            % [WLS Core] Compute weight based on CRLB (Variance = 1/SNR)
            % Sigma^2 = c^2 / (pi^2 * BW^2 * 2 * SNR)
            snr_lin = 10^(snr_db(i)/10);
            sigma_sq = (c^2) / ( (pi^2) * (prs_bw^2) * 2 * snr_lin );
            W(i, i) = 1 / sigma_sq;
        end
        
        % [WLS Core] Solve: (H' * W * H) * dx = H' * W * delta_rho
        delta_x = (H' * W * H) \ (H' * W * delta_rho); 
        
        x_state = x_state + delta_x;
        if norm(delta_x(1:3)) < tol, break; end
    end
    
    est_pos = x_state(1:3)';
    est_clock = x_state(4);
    Q = inv(H' * W * H); % DOP now calculated with weights
    GDOP = sqrt(trace(Q)); 
end