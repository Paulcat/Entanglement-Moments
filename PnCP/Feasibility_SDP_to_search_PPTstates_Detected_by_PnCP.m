% Continue from your initialization
clear all;
d = 3;
load('PnCP_Igor.mat')
load('rho_initialization_dim3');

% Initialize YALMIP variables
yalmip('clear');
rho = sdpvar(d*d, d*d, 'hermitian', 'complex');
assign(rho, delta);  % Starting point from your initialization

% Apply the linear transformation
A = ApplyPnCPSym_left(d, d, rho, phi);

% Instead of using min(eig(A)), we'll use a new variable
t = sdpvar(1,1);  % Variable for the eigenvalue bound

% Define constraints
Constraints = [
    trace(rho) == 1,          % Trace constraint
    rho >= 0,                 % PSD constraint
    A + t*eye(size(A)) >= 0,  % This ensures A has an eigenvalue <= -t
    t >= 0.01                 % We want the negative eigenvalue to be at least -0.01
];

% Simple objective
Objective = 1;

% Set solver options
opts = sdpsettings('solver', 'mosek', 'verbose', 2);

% Solve the problem
sol = optimize(Constraints, Objective, opts);

% Extract results
rho_result = value(rho);
A_result = ApplyPnCPSym_left(d, d, rho_result, phi);
eigenvalues = eig(A_result);

% Display results
fprintf('Optimization completed with status: %s\n', sol.info);
fprintf('Minimum eigenvalue achieved: %f\n', min(real(eigenvalues)));
disp('All eigenvalues:');
disp(eigenvalues);

% Verify constraints
fprintf('Trace of rho: %f\n', trace(rho_result));
fprintf('Minimum eigenvalue of rho: %f\n', min(eig(rho_result)));

