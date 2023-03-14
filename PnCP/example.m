% example of how to use gen_PnCP

n = 3;
m = 3;


% set options
options = struct;
options.mode   = 'gen-many'; % stop searching after finding 'quality' map
options.tol   = 1e-4; % criterion of quality
options.ntest = 50;
options.method = 'klep'; % choice of implementation
options.maxorder = 1; % maximum relaxation order
options.solver = 'mosek';
options.verbose = 0;

% generate PnCP
[phi,delta] = gen_PnCP(n,m,options);
