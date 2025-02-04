clear all;
yalmip('clear');

d = 3;
load('PnCP_Igor.mat');
load('rho_initialization_dim3');  % Load an initial rho solution

% Define YALMIP variable (Hermitian positive semidefinite matrix)
rho = sdpvar(d*d, d*d, 'hermitian', 'complex');
assign(rho, delta);  % Set an initial value for rho

% Define the minimal eigenvalue variable (scalar t)
t = sdpvar(1, 1);

% Define the SDP constraints
F = [rho >= 0, trace(rho) == 1];  % rho should be positive semidefinite and have unit trace

% Constraints to ensure t is the smallest eigenvalue of rho
F = [F, rho - t*eye(d*d) >= 0];  % Ensure that rho - t*I is positive semidefinite

% Define the objective to maximize the smallest eigenvalue t
Objective = -t;  % Maximizing t is equivalent to minimizing the smallest eigenvalue

% Set solver options (using MOSEK, or you can switch to 'sdpt3' or 'sedumi')
options = sdpsettings('solver', 'mosek', 'verbose', 1);

% Solve the optimization problem
sol = optimize(F, Objective, options);

% Display results
if sol.problem == 0
    optimal_rho = value(rho);
    optimal_value = value(Objective);
    
    disp('Optimal solution found:');
    disp('Optimal rho:');
    disp(optimal_rho);
    disp('Optimal value:');
    disp(optimal_value);
else
    disp('Optimization failed. Solver message:');
    disp(sol.info);
end
