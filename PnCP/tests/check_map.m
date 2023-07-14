function [status,obj] = check_map(Phi,n,m,order)
%CHECK_MAP check positivity of a map
%   S = CHECK_MAP(Phi,n,m,order) check if the map Phi: Mn(R) -> Mm(R) is
%   positive by checking positivity of the corresponding polynomial using
%   gloptipoly.
%
%   order specifies the relaxation order
%
%   Flag : 1 : the map is positive,
%          0 : indefinite,
%          -1: the map is not positive

p = map2poly(Phi);

% TODO: can we automatize n, m? --> yes:
mpol('x',n,1);
mpol('y',m,1);

% monomials xi*xj with 1<=i<=j<=n
[A,B] = meshgrid(1:n);
ij = [A(:),B(:)]; A = []; B = [];
%
X = x*x';
X = X(ij(:,2)>=ij(:,1));

% monomials yk*yl with 1<ki<=l<=n
[A,B] = meshgrid(1:m);
kl = [A(:),B(:)]; A = []; B = [];
%
Y = y*y';
Y = Y(kl(:,2)>=kl(:,1));

% declare polynomial
P = X' * p * Y;

Pb = msdp(min(P),order); % TODO: allow to fix relaxation order outside of the function


% code too heavy for local laptop...
prompt = 'Memory required, MATLAB may crash. Do you want to proceed? Y/N [N]: ';
txt = input(prompt,'s');
if isempty(txt)
	txt = 'N';
end
if strcmpi(txt,'Y')
	[status,obj] = msol(Pb);
else
	fprintf('Code did not run\n');
	status = 0;
	obj = -Inf;
end

%if status == 1 && obj > -Inf % collapse + finite min
%	flag = 1;
%elseif status == 1 && obj==-Inf % collapse + no lower bound
%	flag = -1;
%elseif status == 0 && obj > -Inf % no collapse, but finite lower bound
%	flag = 1;
%elseif status == 0 && obj==-Inf % no collapse, no lower bound
%	flag = 0;
%elseif status == -1
%	flag = 0;
%	warning('SDP could not be solved');
%else
%	error('unexpected output from gloptipoly');
%end

end