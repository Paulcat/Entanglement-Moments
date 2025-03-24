function [distance] = Distinguishing_PT_PhiPos_Sets_From_Unitaries(Witness,d)


load('PnCP_Igor.mat')



%% 1st SDP

% Initialize YALMIP variables
yalmip('clear');
rho = sdpvar(d*d, d*d, 'hermitian', 'complex');  

% SDP Constraints
Constraints = [rho >= 0, trace(rho)==1];  % rho is SDP and trace 1

% Ensure the states detected are positive under the action of 1| o phi
Constraints = [Constraints, ApplyPnCPSym_left(d, d, rho, phi)  >= 0];  

% Minimization
Objective = trace(Witness*rho);

% Solver options for MOSEK
options = sdpsettings('solver', 'mosek', 'verbose', 0);

% Solve the problem
optimize(Constraints, Objective, options);

distance_Phiset = value(Objective) ;


%% 2nd SDP

% Initialize YALMIP variables
yalmip('clear');
rho = sdpvar(d*d, d*d, 'hermitian', 'complex');  

% SDP Constraints
Constraints = [rho >= 0, trace(rho)==1];  % rho is SDP and trace 1

% Ensure the states detected are positive under the action of 1| o phi
Constraints = [Constraints, PartialTranspose(rho, 2, 3)  >= 0];  

% Minimization
Objective = trace(Witness*rho);

% Solver options for MOSEK
options = sdpsettings('solver', 'mosek', 'verbose', 0);

% Solve the problem
optimize(Constraints, Objective, options);

distance_PPTset = value(Objective) ;


%%

distance = distance_Phiset -  distance_PPTset ;



end
