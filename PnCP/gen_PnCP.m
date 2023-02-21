function [Maps,Deltas] = gen_PnCP(n,m,options)
%GEN_PNCP Generate PnCP map(s)
%   Successive maps are generate by changing vectors in kernels, not
%   the seed points.
%
%   Supported options:
%    - mode: 
%       'hit-tol':  return one PnCP map satysfing criterion (default)
%       'gen-many': return all maps found, including bad ones
%    - ntest: maximum number of tries for PnCP maps
%    - tolerance
%
%   See also GEN_ONE_PNCP

d = n+m-2;


% process input options % TODO: put in separate function
mode  = getoptions(options,'mode','hit-tol');
ntest = getoptions(options,'ntest',100);
tol   = getoptions(options,'tol',1e-1);

tol_mode = strcmp(mode,'hit-tol'); % check for tolerance mode


% specific initial points and resulting kernels
[~,~,Z]     = Klep_step1_1(n,m,'example');
[~,K,K1]    = Klep_step1_2(n,m,Z);
[~,K_inter] = Klep_step2(n,m,Z);


% dimensions
dK = size(K,2);
d1 = size(K1,2);
di = size(K_inter,2);
K  = reshape(K,[],1,dK); 

Maps    = [];
Deltas  = [];

i = 0;
while i <= ntest
	fprintf('%i ',i);
	if ~mod(i,30)
		fprintf('\n');
	end

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
	[phi,delta]     = gen_one_PnCP(n,m,vf,vh,'verbose',0);
	
	% store results
	Maps(:,:,end+1) = phi;
	Deltas(end+1)   = delta;

	if (delta > tol) && tol_mode
		Maps   = Maps(:,:,end);
		Deltas = Deltas(end);
		break;
	end
	
	i = i+1;
end
fprintf('\n\n');

end

