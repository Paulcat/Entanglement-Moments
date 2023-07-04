function p = map2poly(Phi)
%MAP2POLY returns the polynomial associated with a map
%   P = MAP2POLY(PHI) returns the polynomial
%
%        P(x,y) = y* Phi(xx*) y
%
%   See [Klep et al., 2017] for further details
%
%   The coefficients of the polynomial are returned in the following way: a
%   matrix, lines indexed by (i,j) with 1<=i<=j<=n (x) and columns indexed
%   by (k,l) with 1<=k<=l<=m (y), such that at entry ((i,j),(k,l)) one
%   finds the coefficient of the monomial xi xj yk yl


[m,n] = size(Phi);
n     = sqrt(n);
m     = sqrt(m);

% dimension of (degree 2) bihomogenous polynomial in (n,m) variables
deg = 2;
N = nchoosek(n+deg-1,deg);
M = nchoosek(m+deg-1,deg);

% list of indices (i,j) with 1<=i<=j<=n
[A,B] = meshgrid(1:n);
ij = [A(:),B(:)]; A = []; B = [];
ij(find(ij(:,1)>ij(:,2)),:) = [];
%
[A,B] = meshgrid(1:m);
kl = [A(:),B(:)]; A = []; B = [];
kl(find(kl(:,1)>kl(:,2)),:) = [];

if mod(n,1) || mod(m,1)
	error('Map should be passed as a matrix of dimension (m^2 x n^2)');
end

%H = gensymbasis(n,m);
%nb  = numel(H);

%if nb ~= N*M
%	error('Basis for map has an incorrect number of elements. Please check n and m');
%end

Phi = reshape(Phi,[m,m,n,n]);

p = zeros(N,M);
for a=1:N
	i_j = ij(a,:);
	
	if i_j(1) ~= i_j(2)
		phia = Phi(:,:,i_j(1),i_j(2)) + Phi(:,:,i_j(2),i_j(1));
	else
		phia = Phi(:,:,i_j(1),i_j(2));
	end
	
	for b=1:M
		k_l = kl(b,:);
		
		if k_l(1) ~= k_l(2)
			p(a,b) = phia(k_l(1),k_l(2)) + phia(k_l(2),k_l(1));
		else
			p(a,b) = phia(k_l(1),k_l(2));
		end
	end
end

