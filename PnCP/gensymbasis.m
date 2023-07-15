function [H,Sa,Sb] = gensymbasis(n,m)
%GENSYMBASIS Generate elementary symmetric matrices
%   Detailed explanation goes here

kn = 2; % degree wrt x
km = 2; % degree wrt y

N = nchoosek(n+kn-1,kn); % dimension symmetric / bihomogenous polynomials
M = nchoosek(m+km-1,km);

% list of 1<=i<=j<=n
[I,J] = meshgrid(1:n);
ij = [I(:),J(:)];
ij(ij(:,1)>ij(:,2),:) = [];
%
wa = ones(N,1);
wa(ij(:,1)==ij(:,2)) = 1/2;

% list of 1<=k<=l<=m
[K,L] = meshgrid(1:m);
kl = [K(:),L(:)];
kl(kl(:,1)>kl(:,2),:) = [];
%
wb = 1/2 * ones(N,1);
wb(ij(:,1)==ij(:,2)) = 1;

for a=1:N
	i = ij(a,1);
	j = ij(a,2);
	%
	Ta  = sparse(i,j,wa(a),n,n) + sparse(j,i,wa(a),n,n);
	
	for b=1:M
		k = kl(b,1);
		l = kl(b,2);
		%
		Tb  = sparse(k,l,wb(b),m,m) + sparse(l,k,wb(b),m,m);
		
		Sa{a,b} = Ta;
		Sb{a,b} = Tb;
		H{a,b}  = kron(Ta,Tb);
		
	end
end

end

