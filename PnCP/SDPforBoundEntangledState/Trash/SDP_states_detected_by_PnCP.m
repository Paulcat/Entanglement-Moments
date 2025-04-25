clear all;

%%% WRONG CONSTRAINTS NOT CONVEX


d=3;
load('PnCP_Igor.mat')
load('rho_initialization_dim3');  % Load an initial rho solution
%phi=PnCP_Igor.phi;

% Initialize YALMIP variables
yalmip('clear');
rho = sdpvar(d*d, d*d, 'hermitian', 'complex');  

eps=1e-2;

% Hyperplane that separates states that are positive under 1| o phi from
% other states
W_Ig = 2/3 * eye(9) - delta;


%% SDP Constraints
Constraints = [rho >= 0, trace(rho)==1];  % rho is SDP and trace 1

% Ensure states are detected by phi
Constraints = [Constraints, trace(rho*delta) <= -eps];  

%Constraints = [Constraints, PartialTranspose(rho)  >= 0];  


%%
% Minimization
Objective = 1;

% Solver options for MOSEK
options = sdpsettings('solver', 'mosek', 'verbose', 1);

% Solve the problem
optimize(Constraints, Objective, options);



%% Display results

disp('Best rho found:');
disp(value(rho));
