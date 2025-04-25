% toy example for generating and applying a PnCP

% set path: needs yalmip, mosek, gloptipoly, sedumi
%addpath('../');


% set up problem

% choice of a bipartite decomposition
dA = 4;
dB = 4;

% density matrix
%rho = load('rho3x3.txt');

% dimensions for PnCP map Phi: M_n(C) -> M_m(C)
n = dB; % necessarily n = dB
m = dB; 
% NB: always n = m for us
% but one may also consider m != dB (see [Klep et al.,2017])

% set options (see gen_pncp for details)
options = struct;
options.tol      = 1e-3; % tolerance on PnCP map 'quality'
options.ntest    = 1000; % maximum number of trials
options.method   = 'klep'; % choice of implementation
options.maxorder = 3; % maximum relaxation order
options.solver   = 'mosek';
options.verbose  = 1;

% generate PnCP
nmap  = 1;
PnCPs  = cell(1,nmap);
deltas = zeros(1,20);
for i = 1:nmap
   flag = 0;
   while flag ~= 1
      [phi,delta,flag] = gen_PnCP(n,m,options);
   end

   PnCPs{i}  = phi;
   deltas(i) = delta;
end

% apply PnCP
% MatAfterPnCP = ApplyPnCPSym_right(dA,dB,rho,phi);

%% test positivity

% needs gloptipoly, sedumi
addpath(genpath('~/Workspace/build/gloptipoly3'));
addpath(genpath('/export/home/pcatala/Documents/MATLAB/cvx/')); % for sedumi

% naive, unreliable test
test = test_positive(phi);

% reliable test, using Lasserre's hierarchy
Maps = load('maplist30.mat');
Maps = Maps.map_list;
Map  = Maps{3}; % for example

n = 4;
m = 4;
flag = check_map(Map,n,m,3);








%% obsolete
% small test positivity of polynomial with gloptipoly3

%addpath(genpath('~/Workspace/build/gloptipoly3'));
%addpath(genpath('/export/home/pcatala/Documents/MATLAB/cvx/')); % for sedumi

mpol x1 x2 x3 y1 y2 y3;
P  = ...
	  104*x1^2*y1^2 + 283*x1^2*y2^2 + 18*x1^2*y3^2 - 310*x1^2*y1*y2 ...
	+ 18*x1^2*y1*y3 + 4*x1^2*y2*y3 + 310*x1*x2*y1^2 - 18*x1*x3*y1^2 ...
	- 16*x1*x2*y2^2 + 52*x1*x3*y2^2 + 4*x1*x2*y3^2 - 26*x1*x3*y3^2 ...
	- 610*x1*x2*y1*y2 - 44*x1*x3*y1*y2 + 36*x1*x2*y1*y3 - 200*x1*x3*y1*y3 ...
	- 44*x1*x2*y2*y3 + 322*x1*x3*y2*y3 + 285*x2^2*y1^2 + 16*x3^2*y1^2 ...
	+ 4*x2*x3*y1^2 + 63*x2^2*y2^2 + 9*x3^2*y2^2 + 20*x2*x3*y2^2 ...
	+ 7*x2^2*y3^2 + 125*x3^2*y3^2 - 20*x2*x3*y3^2 + 16*x2^2*y1*y2 ...
	+ 4*x3^2*y1*y2 - 60*x2*x3*y1*y2 + 52*x2^2*y1*y3 + 26*x3^2*y1*y3 ...
	- 330*x2*x3*y1*y3 - 20*x2^2*y2*y3 + 20*x3^2*y2*y3 - 100*x2*x3*y2*y3;

% P    = ...
% 	  104*x(1)^2*y(1)^2 + 283*x(1)^2*y(2)^2 + 18*x(1)^2*y(3)^2 ...
% 	- 310*x(1)^2*y(1)*y(2) + 18*x(1)^2*y(1)*y(3) + 4*x(1)^2*y(2)*y(3) ...
% 	+ 310*x(1)*x(2)*y(1)^2 - 18*x(1)*x(3)*y(1)^2 - 16*x(1)*x(2)*y(2)^2 ...
% 	+ 52*x(1)*x(3)*y(2)^2 + 4*x(1)*x(2)*y(3)^2 - 26*x(1)*x(3)*y(3)^2 ...
% 	- 610*x(1)*x(2)*y(1)*y(2) - 44*x(1)*x(3)*y(1)*y(2) ...
% 	+ 36*x(1)*x(2)*y(1)*y(3) - 200*x(1)*x(3)*y(1)*y(3) ...
% 	- 44*x(1)*x(2)*y(2)*y(3) + 322*x(1)*x(3)*y(2)*y(3) + 285*x(2)^2*y(1)^2 ...
% 	+ 16*x(3)^2*y(1)^2 + 4*x(2)*x(3)*y(1)^2 + 63*x(2)^2*y(2)^2 ...
% 	+ 9*x(3)^2*y(2)^2 + 20*x(2)*x(3)*y(2)^2 + 7*x(2)^2*y(3)^2 ...
% 	+ 125*x(3)^2*y(3)^2 - 20*x(2)*x(3)*y(3)^2 + 16*x(2)^2*y(1)*y(2) ...
% 	+ 4*x(3)^2*y(1)*y(2) - 60*x(2)*x(3)*y(1)*y(2) + 52*x(2)^2*y(1)*y(3) ...
% 	+ 26*x(3)^2*y(1)*y(3) - 330*x(2)*x(3)*y(1)*y(3) - 20*x(2)^2*y(2)*y(3) ...
% 	+ 20*x(3)^2*y(2)*y(3) - 100*x(2)*x(3)*y(2)*y(3);

Pb = msdp(min(P),6);

[status,obj] = msol(Pb);

%% Motzkin polynomial

mpol x y;
P = x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1;

Pb = msdp(min(P),10);
[status,obj] = msol(Pb);
