clear all;

% Define the function to maximize
objectiveFunction = @(alpha) -Distinguishing_PT_PhiPos_Sets(alpha);

Number_of_Optim = 100;

Distance= zeros(1,Number_of_Optim);

parfor i=1:Number_of_Optim
    %disp(i)


%initial guess
alpha0=rand(1,8);
alpha0=alpha0/sum(alpha0);

% Constraints: sum(alpha) = 1 and alpha(i) ≥ 0
Aeq = ones(1, 8);  % Equality constraint: sum(alpha) = 1
beq = 1;
lb = zeros(1, 8);  % Lower bound: alpha(i) ≥ 0
ub = ones(1, 8);   % Upper bound (not strictly needed)

% Optimization options
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');

% Run constrained optimization
[alpha_opt, F_max] = fmincon(objectiveFunction, alpha0, [], [], Aeq, beq, lb, ub, [], options);

% % Display results
% disp('Optimal alpha:');
% disp(alpha_opt);
% disp('Maximum function value:');
% disp(-F_max);  % Convert back to maximization

Distance(i) = -F_max;

end

max(Distance)
