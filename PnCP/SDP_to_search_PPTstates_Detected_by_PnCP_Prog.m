d=3;

phi=PnCP_Igor.phi;

% Initialize YALMIP variables
yalmip('clear');
rho = sdpvar(d*d, d*d, 'hermitian', 'complex');  % Hermitian variable

% Apply the custom operation
A = ApplyPnCPSym_left(d, d, rho, phi); % (1| o Phi) (rho) 

% Define eigenvalue minimization
% Introduce a scalar variable for the minimal eigenvalue
t = sdpvar(1, 1);

% SDP Constraints
Constraints = [rho >= 0, trace(rho) == 1, A - t*eye(size(A)) >= 0]; %(1| o Phi) (rho)  - t Id > 0 

% Objective: Minimize the smallest eigenvalue
Objective = 1/t;

% Solver options for MOSEK
options = sdpsettings('solver', 'mosek', 'verbose', 1);

% Solve the problem
optimize(Constraints, Objective, options);

% Display results
optimal_rho = value(rho);
optimal_value = value(Objective);

disp('Optimal rho:');
disp(optimal_rho);
disp('Optimal value:');
disp(optimal_value);