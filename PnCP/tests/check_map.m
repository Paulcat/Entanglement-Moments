function flag = check_map(Phi,n,m,max_order)
%CHECK_MAP check positivity of a map
%   S = CHECK_MAP(Phi,n,m,order) checks if the map Phi: Sn -> Sm is PnCP by
%   checking nonnegatvity and sos-ity of the corresponding polynomial. 
%
%   Inputs:
%      - n,m:       input and output dimensions of the map (integers)
%      - Phi:       map (m^2 x n^2 matrix)
%      - max_order: maximum relaxation order for the gloptipoly test
%
%   Flag : 1 : the map is PnCP,
%          0 : indefinite,
%          -1: the map is not PnCP


p = map2poly(Phi);

% helpers
[A,B] = meshgrid(1:n);
[C,D] = meshgrid(1:m);
%
ij = [A(:),B(:)]; % list of (i,j)
kl = [C(:),D(:)]; % list of (k.l)
%
A = []; B = []; C = []; D = [];

% 1- checking sos positivity of p (Yalmip)
% ***
x = sdpvar(n,1);
y = sdpvar(m,1);

% monomials xi*xj with i<=j
X = x*x';
X = X(ij(:,2)>=ij(:,1));

% monomials yk*yl with k<=l
Y = y*y';
Y = Y(kl(:,2)>=kl(:,1));

% declare bihomogeneous polynomial for Yalmip
P1 = X' * p * Y;

% check sos
[sol,res] = is_sos(P1,'solver','mosek');




% 2 - checking lower bound on polynomial (gloptipoly)
% ***
clear x X y Y
mpol('x',n,1);
mpol('y',m,1);


% monomials, as above
X = x*x';
X = X(ij(:,2)>=ij(:,1));
%
Y = y*y';
Y = Y(kl(:,2)>=kl(:,1));

% declare polynomial for gloptipoly
P = X' * p * Y;


status = -1;
order  = 2; % at least half of the polynomial degree

while status ~= 1 && order <= max_order
	
	% declare problem
	Pb = msdp(min(P),order);
	
	% code may be heavy for local laptop...
	fprintf('Relaxation order: %i\n', order);
	%prompt = 'MATLAB may crash for large orders. Do you want to proceed? Y/N [N]: ';
	%txt = input(prompt,'s');
	%if isempty(txt)
	%	txt = 'N';
	%end
	
	%if strcmpi(txt,'Y')
		[status,obj] = msol(Pb);
	%else
	%	fprintf('Code did not run\n');
	%	status = 0;
	%	obj = -Inf;
	%end
	
	order = order+1;
end

% summary log
fprintf('\n\n *** RESULT ***\n');
if sol.problem == 0
	% sos positivity is satisfied (up to tolerance)
	fprintf('Map is completely positive\n');
	fprintf('\t residual = %d\n', res);
	
	flag = 2;
	
elseif status ~= -1 && obj > -Inf
	% sos is not satisfied
	fprintf('%20s\n', 'Map is not completely positive');
	fprintf('%+20s : %s\n', 'info (Y)', sol.info);
	
	% nonnegativity is satisfied
	fprintf('%20s\n', 'Map is positive');
	fprintf('%+20s : %i (0: no collapse, 1: collapse)\n', 'status (G)', status);
	fprintf('%+20s : %d\n', 'min/bound (G)', obj);
	
	flag = 1;
	
elseif status == 1 && obj == -Inf
	% sos is not satisfied
	fprintf('%20s\n', 'Map is not completely positive'); 
	fprintf('%+20s : %s\n', 'info (Y)', sol.info);
	
	% nonnegativity is not satisfied
	fprintf('%20s\n', 'Map is not positive');
	fprintf('%+20s : %i (1: collapse)\n', 'status (G)', status);
	fprintf('%+20s : %d\n', 'min (G)', obj);

   flag = -1;
	
else
	% sos is not satisfied
	fprintf('%20s\n', 'Map is not completely positive');
	fprintf('%+20s : %s\n', 'info (Y)', sol.info);
	
	% we cannot conclude on nonnegativity
	fprintf('%20s\n', 'Nonnegativity could not be checked');
	fprintf('%+20s : %i (0: no collapse, -1: not solved)\n', 'status (G)', status);
	fprintf('%+20s : %d\n', 'min/bound (G)', obj);

   flag = 0;
end
fprintf('\n');
