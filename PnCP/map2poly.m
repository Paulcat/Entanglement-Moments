function p = map2poly(Phi)
%MAP2POLY returns the polynomial associated with a map
%   P = MAP2POLY(PHI) returns the polynomial
%
%        P(x,y) = y'PHI(xx')y
%
%   associated with the map PHI (see [Klep et al., 2017] for details). The
%   input PHI is a matrix of size m^2 x n^2.
%
%   P is a matrix indexed by r=(i,j) (1<=i<=j<=n) and s=(k,l) (1<=k<=l<=m),
%   in colexicographic order, such that P(r,s) is the coefficient of the
%   monomial xi*xj*yk*yl.
%
%   See also poly2map

[m,n] = size(Phi);
n     = sqrt(n);
m     = sqrt(m);


if mod(n,1) || mod(m,1)
	error('Map should be passed as a matrix of dimension (m^2 x n^2)');
end

% dimension of (2,2)-bihomogenous polynomials in (n,m) variables
deg = 2;
N = nchoosek(n+deg-1,deg); % N = n(n+1)/2
M = nchoosek(m+deg-1,deg);

% bases for symmetric matrices
Sn = gensymbasis(n);
Sm = gensymbasis(m);

p = zeros(N,M);
for a=1:N
	phia = reshape(Phi * Sn{a}(:),[n,n]);
	
	for b=1:M
		p(a,b) = trace(Sm{b}'*phia);
	end
	
end


% obsolete

% list of indices (i,j) with 1<=i<=j<=n
% [A,B] = meshgrid(1:n);
% ij = [A(:),B(:)]; A = []; B = [];
% ij(find(ij(:,1)>ij(:,2)),:) = []; % bijection from Delta to D (see notes)
% 
% [A,B] = meshgrid(1:m);
% kl = [A(:),B(:)]; A = []; B = [];
% kl(find(kl(:,1)>kl(:,2)),:) = [];
% 
% Phi2 = reshape(Phi,m,m,n,n);


end

