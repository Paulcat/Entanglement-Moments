function [H,Hn] = gen_sym_basis(n,m)
%GEN_SYM_BASIS Generate elementary symmetric matrices
%   GEN_SYM_BASIS(N,M) returns the list of matrices units H(i,j,k,l) such
%   that
%
%        v(x,y) * v(x,y)' =    Sum      H((i,j),(k,l)) x(i)x(j) y(k)y(l)
%                           i<=j, k<=l
%
%   where v is the vector of (1,1)-bihomegeneous monomials
%
%        v(x,y) = ( x . y )
%
%   v is of size (n*m), and H is a tensor of size (n*m x n*m x N x M)
%
%   The system {H(i,j,k,l)} is a basis for (PnCP) maps (seen as tensors)
%   (which are in bijection with the set of (2,2)-bihomogeneous
%   polynomials)


kn = 2;
km = 2;

N = nchoosek(n+kn-1,kn); % dimension of (kn)-homogeneous polynomials
M = nchoosek(m+km-1,km);

% to put entries of (PnCP) map seen as n*m x n*m matrix in correlation with
% coefficients of polynomial, enumerate block-wise upper triangle from left
% to right, bottom to top. Doing so, the "diagonal" elements can be
% assigned numbers as follow
dn = cumsum([1,n:-1:2]);
dm = cumsum([1,m:-1:2]);

% assign numbers to remaining off-diagonal entries
on = setdiff(1:N,dn);
om = setdiff(1:M,dm);

% TODO: a specific ordering is assumed. Make it more flexible?
% TODO: unnecessary complicated code due to this choice of ordering?
% (require weird procedure to identify diagonal elements...)


%H = cell(N,M);
H  = zeros(n*m,n*m,N,M);
Hn = zeros(n*m,n*m,N,M);

for i=1:N
	% separate between diagonal/off-diagonal elements
	idO = find(on==i);
	idD = find(dn==i);
	
	if idO
		% elementary extraction matrix
		theta_i          = sym_off(n,idO); % Eij
	else
		theta_i          = zeros(n);
		theta_i(idD,idD) = 1; % Eii
	end
	
	for j=1:M
		% same thing
		idO = find(om==j);
		idD = find(dm==j);
		
		if idO
			% elementary matrix
			theta_j          = sym_off(m,idO);
		else
			theta_j          = zeros(m);
			theta_j(idD,idD) = 1;
		end
		
		Hij = kron(theta_i,theta_j);
		H (:,:,i,j) = Hij;
		Hn(:,:,i,j) = Hij / norm(Hij,'fro')^2;
		
		%H{i,j} = kron(theta_i,theta_j);
		%C(j,i)  = trace(thetaij'*v);
	end
end


end

function Ho = sym_off(n,s)
%HANKEL_OFF Off-diagonal symmetric elementary matrix of size n
%   Detailed explanation goes here

offdiag = ones(n);
offdiag(logical(eye(n))) = 0;
[A,B] = meshgrid(0:n-1);

Ho = double(A+B==s) .* offdiag;
end



