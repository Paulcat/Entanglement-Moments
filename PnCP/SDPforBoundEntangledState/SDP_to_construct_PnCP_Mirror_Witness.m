clear all;

d=3;
load('PnCP_Igor.mat')
load('rho_initialization_dim3');  % Load an initial rho solution
%phi=PnCP_Igor.phi;

% Initialize YALMIP variables
yalmip('clear');
rho = sdpvar(d*d, d*d, 'hermitian', 'complex');  

W_Ig = 2/3 * eye(9) - delta;

W_nonopt = 1 * eye(9) - delta;

W_MirrorNonOpt = delta;


% Apply the custom operation
A = ApplyPnCPSym_left(d, d, rho, phi); % Compute A = phi(rho)


%% SDP Constraints
Constraints = [rho >= 0, trace(rho)==1];  % rho is SDP and trace 1

% Ensure the states detected are positive under the action of 1| o phi
Constraints = [Constraints, A  >= 0];  


%%
% Minimization
Objective = 1-trace(rho*delta);

% Solver options for MOSEK
options = sdpsettings('solver', 'mosek', 'verbose', 1);

% Solve the problem
optimize(Constraints, Objective, options);



%% Display results

disp('Objective function');
disp(value(Objective));

% disp('Best rho found:');
% disp(value(rho));
