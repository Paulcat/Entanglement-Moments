function S = gensymbasis(n)
%GENSYMBASIS Generate elementary symmetric matrices
%   B = GENSYMBASIS(N) returns the list of elementary symmetric matrices
%
%                  | E(i,j) + E(j,i) if i<j
%        S(i,j) = <
%                  | E(i,i) if i=j
%
%    The matrices are returned in colexicographical order, i.e.
%    B{1} = S(1,1), B{2} = S(1,2), ..., B{n(n+1)/2} = S(n,n).
%
%    Internal use only:
%    See also MAP2POLY

k = 2; % degree wrt x

N = nchoosek(n+k-1,k); % dimension symmetric / bihomogenous polynomials

% list of 1<=i<=j<=n
[I,J] = meshgrid(1:n);
ij    = [I(:),J(:)];
ij(ij(:,1)>ij(:,2),:) = [];
%
wa = ones(N,1);
wa(ij(:,1)==ij(:,2)) = 1/2;

S = cell(N,1);
for a=1:N
	i = ij(a,1);
	j = ij(a,2);
	
	S{a} = sparse(i,j,wa(a),n,n) + sparse(j,i,wa(a),n,n);
end

end

