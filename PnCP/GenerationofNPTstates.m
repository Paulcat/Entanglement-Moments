% Define the SDP variable (a 9x9 matrix)
X = sdpvar(9, 9, 'hermitian');

% Objective: maximize 1 (trivial objective)
objective = 1;

% Constraint: X must be positive semidefinite
constraints = [X >= 0, trace(X) == 1];

% Set options for the solver (MOSEK)
ops = sdpsettings('solver', 'mosek');

% Solve the optimization problem
sol = optimize(constraints, objective, ops);  

% Display the results
if sol.problem == 0
    disp('Optimization successful!');
    disp('Optimal value for the objective:');
    disp(value(objective));
    
    % Retrieve the optimal value of X
    optimal_X = value(X);
    disp('Optimal value of X:');
    disp(optimal_X);
    
else
    disp('Optimization failed.');
    disp('Solver message:');
    disp(sol.info);
end
