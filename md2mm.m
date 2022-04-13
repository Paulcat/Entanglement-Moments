function M = md2mm(D,m,n,E)
%MD2MM Compute moment matrix from density matrix of a bipartite state
%   M = MD2MM(D,M,N,MON) returns the matrix such that
%
%		M[k,l] = Tr{ f_k' * f_l * D }
%
%	 where 
%
%		f_i = {Adag^E(i,1) * A^E(i,2)} x {Bdag^E(i,3) * B^E(i,4)},
%
%	 E is a L x 4 matrix whose lines contain the quadruplet exponents for 
%	 (Adga, A, Bdga, B), the respective creation/destruction operators.
%
%	 M and N specify the dimension of the two Hilbert spaces.
%
% 	 Requires MULTIPROD, Paolo de Leva (2022) MATLAB Central File Exchange

L = size(E,1); % size of output moment matrix
q = m*n;

[a,ad] = quantop(m);
[b,bd] = quantop(n);

C 	 = num2cell(E);
Phi = zeros(q,q,L);
%
for t=1:L %TODO: maybe doable without for loop, using *tensorprod* introduced in R2022a...
	[i,j,k,l]  = deal(C{t,:}); % exponents of monomials involved in moment matrix
	Phi(:,:,t) = kron(ad^i * a^j, bd^k * b^l);
end

PhiS = conj(permute(Phi,[2 1 3])); % Hermitian transpose of operators
Phi  = reshape(Phi,[q,q,1,L]);


X = multiprod( multiprod(PhiS,Phi), D ); % matrix products along the first 2 dimensions
%TODO: is that better than for loop?
% if >R2020b, multiprod may be replaced with *pagemtimes*
X = reshape(X,[q,q,L^2]);

M = zeros(L^2,1);
for i=1:L^2
	M(i) = trace(X(:,:,i));
end
M = reshape(M,[L,L]);

% alternative, probably worse
%M = arrayfun(@(i)trace(X(:,:,i),1:L^2);
%M = reshape(M,[L,L]);


end

