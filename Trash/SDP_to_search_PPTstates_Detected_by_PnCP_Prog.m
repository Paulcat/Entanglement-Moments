clear all;

d=3;
load('PnCP_Igor.mat')
load('rho_initialization_dim3');  % Load an initial rho solution
%phi=PnCP_Igor.phi;

% Initialize YALMIP variables
yalmip('clear');
rho = sdpvar(d*d, d*d, 'hermitian', 'complex');  
assign(rho, delta);  % Starting point of the optimisation


% Apply the custom operation
A = ApplyPnCPSym_left(d, d, rho, phi); % Compute A = phi(rho)

% Define the minimal eigenvalue variable (scalar t)
t = sdpvar(1, 1);

%% SDP Constraints
Constraints = [rho >= 0, 0.99 <= trace(rho) <= 1.01];  % rho is SDP and trace 1

% Ensure `t` is the smallest eigenvalue of A
Constraints = [Constraints, A - t*eye(size(A)) >= 0];  

% Ensure that `t < 0` to force A to have a negative eigenvalue
Constraints = [Constraints, t <= -1];  

%Constraints = [Constraints, rho >= 1e-6 * eye(size(rho))]; % so that rho is not a 0 state

%%
% Objective: Minimize the smallest eigenvalue
Objective = t;

% Solver options for MOSEK
options = sdpsettings('solver', 'mosek', 'verbose', 1);

% Solve the problem
optimize(Constraints, Objective, options);

% Renormalization
rho_optimal = value(rho);
%rho_optimal=rho_optimal-min(eig((rho_optimal)))*eye(size(rho_optimal));
rho_optimal = rho_optimal / trace(rho_optimal);  % Rescale so trace = 1


%% Display results

%optimal_value = value(Objective);

disp(['Optimal t: ', num2str(value(t))]);

disp('rho optimal eigenvalues:');
disp(eig(rho_optimal));

disp('Eigenvalues of optimal rho after application of PnCP map:');
temp_map=ApplyPnCPSym_left(d, d, rho_optimal, phi);
disp(eig(temp_map));