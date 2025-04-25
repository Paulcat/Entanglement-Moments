function [best_distance,best_alpha] = Optimize_GellMann_Coefficients(rho, num_trials)
    d = 3;
    num_coeffs = 2*d * d; % 9 coefficients

    % Set up optimization options
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off', ...
                           'MaxIterations', 1000, 'TolFun', 1e-10);

    % Preallocate arrays to store results
    distances = zeros(num_trials, 1);
    alpha_matrix = zeros(num_trials, num_coeffs);

    % Run the optimization in parallel
    parfor trial = 1:num_trials
        % Random initial guess
        alpha0 = randn(num_coeffs, 1) * 0.1; 
        
        % Perform optimization
        [alpha_opt, distance] = fminunc(@(alpha) Finding_coefficicents_from_Gellmann_Basis(alpha, rho), alpha0, options);

        % Store results
        distances(trial) = distance;
        alpha_matrix(trial, :) = alpha_opt';
    end

    % Find the best result
    [best_distance, best_index] = min(distances);
    best_alpha = alpha_matrix(best_index, :)';
end

