function [phi,delta,flag] = gen_PnCP(n,m,options,varargin)
%GEN_PNCP Generate a PnCP map following [Klep, 2017]
%   PHI = GEN_PNCP(n,m,options) returns a (PnCP) map.
%
%   The map PHI is a matrix of size (m*m x n*n). It can be applied to a
%   matrix A of size (nxn) by doing PHI * A(:), and reshaping.
%
%   PHI is found by generating maps (successive random initialization)
%   until a criterion is met (sufficient deviation from the sos property).
%
%   Supported options:
%    - ntest (100)      : maximum number of tries (re-initializations)
%    - tol (1e-2)       : tolerance for PnCP map (on deviation from sos)
%    - solver (mosek)   : solver for the SDP program
%    - toolbox (yalmip) : yalmip or (todo: cvx)
%    - verbose (1)
%
%   PHI = GEN_PNCP(N,M,OPTIONS,'EXAMPLE') simply runs the example 4.6 in
%   [Klep et al., 2017].
%
%   See also GEN_ONE_MAP


% process input options % TODO: put in separate function?
nmax    = getoptions(options,'ntest',100);
maxor   = getoptions(options,'maxorder',2);
tol_de  = getoptions(options,'tol',1e-2);
method  = getoptions(options,'method','klep');
solver  = getoptions(options,'solver','mosek');
toolbox = getoptions(options,'toolbox','yalmip');
verbose = getoptions(options,'verbose',1);



% PhiCells={};
% 
% parfor i=1:nmax
i= 0;
while i<=nmax
	%Klep setps 1 and 2
	if nargin < 4
		Z  = Klep_step1_1(n,m);
		vh = Klep_step1_2(n,m,Z);
		vf = Klep_step2(n,m,Z);
	else
		Z  = Klep_step1_1(n,m,'example');
		vh = Klep_step1_2(n,m,Z,'example');
		vf = Klep_step2(n,m,Z,'example');
	end

    Z  = Klep_step1_1(n,m);
	vh = Klep_step1_2(n,m,Z);
	vf = Klep_step2(n,m,Z);

	% Klep step 3: compute PnCP map
	[phi,delta,info]  = gen_one_map(n,m,vf,vh,...
		'verbose',verbose,...
		'maxorder',maxor,...
		'solver',solver,...
		'toolbox',toolbox,...
		'method',method,...
		'tolerance',tol_de);
	
	if nargin >= 4
		vfr = reshape(vf,n*m,n*m);
		vhr = reshape(vh,n*m,[]);
		delta = 2;
		phi = delta*vfr + (vhr*vhr');
		clear vfr vhr
		phi = reshape(phi,[m,n,m,n]);
		phi = permute(phi,[1,3,2,4]);
		phi = reshape(phi,[m*m,n*n]);
	end

% 	if info.success
%         disp(i)
% 		PhiCells{i}=phi;
%         fprintf('\t map found: residual=%d, delta=%d, sdpflag = %i\n',...
% 		info.res,delta,info.flag_sol)
%     else
%         fprintf('\t no satisfying map was found, returning 0\n');
% 	end
	
    i = i+1;
end
fprintf('\n\n');

if i==nmax+1
	% no satisfying map was found
	fprintf('\t no satisfying map was found, returning 0\n');
	phi = 0;
else
	fprintf('\t map found: residual=%d, delta=%d, sdpflag = %i\n',...
		info.res,delta,info.flag_sol)
end

additional post-check
flag = check_map(phi,n,m,3);


end

