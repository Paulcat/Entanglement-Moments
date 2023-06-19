% example of how to use gen_PnCP

n = 2;
m = 4;


% set options
options = struct;
options.mode   = 'hit-tol'; % stop searching after finding 'quality' map
options.tol   = 1e-4; % criterion of quality
options.ntest = 50;
options.method = 'klep'; % choice of implementation
options.maxorder = 1; % maximum relaxation order
options.solver = 'mosek';
options.verbose = 1;

% generate PnCP
[phi,delta] = gen_PnCP(n,m,options);
