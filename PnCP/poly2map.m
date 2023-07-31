function Phi = poly2map(p)
%POLY2MAP returns the map associated with a polynomial
%   Phi = POLY2MAP(p) returns the map such that
%
%          p(x,y) = y'Phi(xx')y
%
%   associated with the map Phi (see [Klep et al., 2017] for details). The
%   input p is a matrix of size n(n+1)/2 x m(m+1)/2.
%
%   See also map2poly


[N,M] = size(p);
n     = (-1+sqrt(1+8*N))/2;
m     = (-1+sqrt(1+8*M))/2;


warning('poly2map might not work if n~=m. Here n=%i, m=%i',n,m);

% indices i<=j
[A,B] = meshgrid(1:n);
ij = [A(:),B(:)]; A = []; B = [];
ij = ij(ij(:,2)>=ij(:,1),:);
%
[A,B] = meshgrid(1:m);
kl = [A(:),B(:)]; A = []; B = [];
kl = kl(kl(:,2)>=kl(:,1),:);


% bases for symmetric matrices
Sn = gensymbasis(n);
Sm = gensymbasis(m);

Phi = zeros(n*m,n*m);
for a=1:N
	% normalize
	Sa = Sn{a};
	Sa = Sa / norm(Sa,'fro')^2;
	
	for b=1:M
		% normalize
		Sb = Sm{b};
		Sb = Sb / norm(Sb,'fro')^2;
		
		Sab = kron(Sa,Sb);
		
		Phi = Phi + p(a,b) * Sab;
	end
end

% put in correct shape
Phi = reshape(Phi,[m,n,m,n]);
Phi = permute(Phi,[1,3,2,4]);
Phi = reshape(Phi,[m^2,n^2]);



end

