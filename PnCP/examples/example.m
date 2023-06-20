% toy example for generating and applying a PnCP

% set path: needs yalmip and mosek
%addpath('../');
%addpath(genpath('/export/home/pcatala/Documents/MATLAB/cvx/mosek'));
%addpath(genpath('~/Workspace/build/yalmip'))

% set up problem
% rho = load('rho3x3.txt');
rho = rand(8);

% choice of a bipartite decomposition
dA = 2;
dB = 4;

% dimensions for PnCP map Phi: M_n(C) -> M_m(C)
n = dA; % necessarily n = dB
m = dB; 
% NB: always n = m for us
% but one may also consider m != dB (see [Klep et al.,2017])

% set options (see gen_pncp for details)
options = struct;
options.mode   = 'hit-tol'; % 'hit-tol' or 'gen_many'
options.tol   = 1e-4; % tolerance on PnCP map 'quality'
options.ntest = 50; % maximum number of trials
options.method = 'klep'; % choice of implementation
options.maxorder = 1; % maximum relaxation order
options.solver = 'mosek';
options.verbose = 1;

% generate PnCP
[phi,delta] = gen_PnCP(n,m,options);

% apply PnCP
MatAfterTest = ApplyPnCP(dA,dB,rho,phi);
