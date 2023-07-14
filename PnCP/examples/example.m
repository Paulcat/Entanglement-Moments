% toy example for generating and applying a PnCP

% set path: needs yalmip and mosek
%addpath('../');
%addpath(genpath('/export/home/pcatala/Documents/MATLAB/cvx/mosek'));
%addpath(genpath('~/Workspace/build/yalmip'))

% set up problem

% choice of a bipartite decomposition
dA = 3;
dB = 3;

% density matrix
% rho = load('rho3x3.txt');
rho = rand(dA*dB);

% dimensions for PnCP map Phi: M_n(C) -> M_m(C)
n = dB; % necessarily n = dB
m = dB; 
% NB: always n = m for us
% but one may also consider m != dB (see [Klep et al.,2017])

% set options (see gen_pncp for details)
options = struct;
options.tol      = 1e-3; % tolerance on PnCP map 'quality'
options.ntest    = 50; % maximum number of trials
options.method   = 'klep'; % choice of implementation
options.maxorder = 1; % maximum relaxation order
options.solver   = 'mosek';
options.verbose  = 1;

% generate PnCP
[phi,delta] = gen_PnCP(n,m,options);

% apply PnCP
MatAfterPnCP = ApplyPnCP_right(dA,dB,rho,phi);

% test positivity
test = test_positive(phi);


