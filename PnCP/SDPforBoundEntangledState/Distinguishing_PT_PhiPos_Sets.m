function [distance] = Distinguishing_PT_PhiPos_Sets(alpha)

%% Define the more general witness
% Define the 8 Gell-Mann matrices
d=3;
GellMann = zeros(d, d, d*d);

GellMann(:,:,1) = [0, 1, 0; 1, 0, 0; 0, 0, 0];                % lambda1
GellMann(:,:,2) = [0, -1i, 0; 1i, 0, 0; 0, 0, 0];             % lambda2
GellMann(:,:,3) = [1, 0, 0; 0, -1, 0; 0, 0, 0];               % lambda3
GellMann(:,:,4) = [0, 0, 1; 0, 0, 0; 1, 0, 0];                % lambda4
GellMann(:,:,5) = [0, 0, -1i; 0, 0, 0; 1i, 0, 0];             % lambda5
GellMann(:,:,6) = [0, 0, 0; 0, 0, 1; 0, 1, 0];                % lambda6
GellMann(:,:,7) = [0, 0, 0; 0, 0, -1i; 0, 1i, 0];             % lambda7
GellMann(:,:,8) = (1/sqrt(3)) * [1, 0, 0; 0, 1, 0; 0, 0, -2]; % lambda8
GellMann(:,:,9) = eye(3);                                     % 1|

% Build the linear combination
Unitary=zeros(3);
for i = 1:d*d
    Unitary = Unitary + alpha(i) * GellMann(:,:,i);
end

Witness=kron(Unitary,Unitary);

%%

load('PnCP_Igor.mat')
%phi=PnCP_Igor.phi;




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
