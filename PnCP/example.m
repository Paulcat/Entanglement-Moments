% example for generating and applying a PnCP


% set up problem
% rho = load('rho3x3.txt');
rho = rand(8);

% choice of a bipartite decomposition
dA = 2;
dB = 4;


% we always have n=m in our problems (cf [Klep,2017] for details on n and m)
n = dB; % necessarily n = dB
m = dB; % from [Klep,2017] one may consider m != dB (but not us)


% set options
options = struct;
options.mode   = 'hit-tol';	% 'hit-tol' : stop searching after finding 'quality' map
										% 'gen-many': generate many maps, whithout guarantee of quality
options.tol   = 1e-4; % criterion on delta (PnCP <-> sos + delta * perturbation)
options.ntest = 50;
options.method = 'klep'; % choice of implementation
options.maxorder = 1; % maximum relaxation order
options.solver = 'mosek';
options.verbose = 1;

% generate PnCP
[phi,delta] = gen_PnCP(n,m,options);

% apply PnCP
MatAfterTest = ApplyPnCP(dA,dB,rho,phi);
