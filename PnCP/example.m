% example of how to use gen_PnCP

n = 3;
m = 3;


% set options
options = struct;
options.mode   = 'hit-tol'; % stop searvhing after finding 'quality' map
options.tol   = 1e-2; % criterion of quality
options.ntest = 50;

% generate PnCP
[phi,delta] = gen_PnCP(n,m,options);
