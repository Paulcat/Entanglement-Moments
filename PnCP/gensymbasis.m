function [H,Ea,Eb,Hn] = gensymbasis(n,m)
%GENSYMBASIS Generate elementary symmetric matrices
%   H = GENSYMBASIS(n,m) returns the list of matrices units H((i,j),(k,l))
%   such that
%
%        v(x,y) * v(x,y)' =    Sum      H((i,j),(k,l)) xi*xj * yk*yl
%                           i<=j, k<=l
%
%   where v(x,y) := (x1*y1, x1*y2, ..., xn*ym). The matrices H((i,j),(k,l))
%   form a basis for maps (which are in bijection with (2,2)-bihomogeneous
%   polynomials, see [Klep et al., 2017]).
%
%   Colexicographical ordering is assumed on the index r=(i,j) and s=(k,l).
%   H is a cell array, H{r,s} is of size n*m x n*m.
%
%   [H,SA,SB] = GENSYMBASIS(N,M) also returns the symmetric matrices
%   SA{r,s} (nxn) and SB{r,s} (mxm) such that H{r,s} = SA{r,s} x SB{r,s}.
%   SA(:,1) (resp. SB(1,:)) gives the canonical basis for symmetric
%   matrices of size n (resp. m).
%
%   [H,SA,SB,Hn] = GENSYMBASIS(N,M) also returns the normalized basis
%   Hn{a,b} = H{a,b} / norm(H{a,b}).


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


H  = cell(N,M);
Ea = cell(N,M);
Eb = cell(N,M);
%H  = zeros(n*m,n*m,N,M);
Hn = zeros(n*m,n*m,N,M);

for i=1:N
	% separate between diagonal/off-diagonal elements
	idO = find(on==i);
	idD = find(dn==i);
	
	if idO
		% elementary extraction matrix
		theta_i          = sym_off(n,idO); % Eij
	else
		%theta_i          = zeros(n);
		%theta_i(idD,idD) = 1; % Eii
		theta_i = sparse(idD,idD,1,n,n);
	end
	
	for j=1:M
		% same thing
		idO = find(om==j);
		idD = find(dm==j);
		
		if idO
			% elementary matrix
			theta_j          = sym_off(m,idO);
		else
			%theta_j          = zeros(m);
			%theta_j(idD,idD) = 1;
			theta_j = sparse(idD,idD,1,m,m);
		end
		
		Hij = kron(theta_i,theta_j);
		
		% store results
		%H (:,:,i,j) = Hij;
		H{i,j} = Hij;
		Ea{i,j} = theta_i;
		Eb{i,j} = theta_j;
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

Ho = sparse(double(A+B==s) .* offdiag);
end



