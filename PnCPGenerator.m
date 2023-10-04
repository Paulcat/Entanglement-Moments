% Generate a list of PnCPs

% bipartite dimensions
dA = 4;
dB = 2;

% dimensions for PnCP map Phi: M_n(C) -> M_m(C)
n = dB;
m = dB;

% set options (see gen_pncp for details)
options = struct;
options.tol      = 1e-3; % tolerance on PnCP map 'quality'
options.ntest    = 10; % maximum number of trials
options.method   = 'klep'; % choice of implementation
options.maxorder = 1; % maximum relaxation order
options.solver   = 'mosek';
options.verbose  = 0;

% generate and store PnCPs
nmaps = 20; % number of maps
PnCPs = {};
for i=1:nmaps
   flag  = 0;
   while flag ~= 1
      [phi,delta,flag] = gen_PnCP(n,m,options);
   end
   PnCPs{i} = phi;
end