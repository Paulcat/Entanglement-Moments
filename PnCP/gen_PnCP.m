function [Maps,Deltas] = gen_PnCP(n,m,options)
%GEN_PNCP Generate PnCP map(s)
%   Successive maps are generate by changing initial seed points, as
%   well as vectors in kernels.
%
%   Supported options:
%    - mode: 
%       'hit-tol' (default): return one PnCP map satysfing criterion
%       'gen-many'         : return all maps found including bad ones
%    - ntest (100)         : maximum number of tries for PnCP maps
%    - tol (1e-2)          : tolerance for the hit-tol mode
%    - solver (sdpt3)      : solver for the SDP program
%    - toolbox (yalmip)    : yalmip or (todo: cvx)
%    - verbose (1)
%
%   See also GEN_ONE_PNCP

d = n+m-2;


% process input options % TODO: put in separate function
mode    = getoptions(options,'mode','hit-tol');
nmax    = getoptions(options,'ntest',100);
maxor   = getoptions(options,'maxorder',2);
tol     = getoptions(options,'tol',1e-1);
method  = getoptions(options,'method','klep');
solver  = getoptions(options,'solver','sdpt3');
toolbox = getoptions(options,'toolbox','yalmip');
verbose = getoptions(options,'verbose',1);

tol_mode = strcmp(mode,'hit-tol'); % check for tolerance mode


Maps   = [];
Deltas = [];
i      = 0;
while i<=nmax
	%fprintf('%i ',i);
	%if ~mod(i+1,30)
	%	fprintf('\n');
	%end

	% random initial points and resulting kernels
	% [~,~,Z]     = Klep_step1_1(n,m,'example'); % specific points
	[~,~,Z]     = Klep_step1_1(n,m);
	[~,K,K1]    = Klep_step1_2(n,m,Z);
	[~,K_inter] = Klep_step2(n,m,Z);

	% dimensions
	dK = size(K,2);
	d1 = size(K1,2);
	di = size(K_inter,2);
	K  = reshape(K,[],1,dK);

	% d random linear combination in K
	al = rand(1,d,dK);
	vj = sum(al.*K,3); % % each column defines the linear form <vj',z>
	vj = reshape(vj,n,m,d); % row and columns are switched (cf Klep_step2.m)
	% (reshape is probably useless)

	% one in K1
	be = rand(1,d1);
	v0 = sum(be.*K1,2);
	v0 = reshape(v0,m,n); % rows and columns switched

	% all vectors
	vh = cat(3,v0,vj);

	% one in K_inter
	ga = rand(1,di);
	vf = sum(ga.*K_inter,2);
	vf = reshape(vf,[m,n,m,n]); % TODO: check if dimensions are in correct order!

	% compute PnCP map
	[phi,delta,~,flag]  = gen_one_PnCP(n,m,vf,vh,...
		'verbose',verbose,'maxorder',maxor,'solver',solver,'toolbox',toolbox,...
		'method',method,'tolerance',tol);
	
	% store results
	Maps(:,:,end+1) = phi;
	Deltas(end+1)   = delta;

	if flag && tol_mode
		Maps   = Maps(:,:,end);
		Deltas = Deltas(end);
		break;
	end
	
	i = i+1;
end
%fprintf('\n\n');

if i==nmax+1 && tol_mode
	% no satisfying map was found
	warning('no satisfying map was found, returning random map');
	Maps = rand(m^2,n^2);
end

end

