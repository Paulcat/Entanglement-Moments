function [v,K_inter] = Klep_step2(n,m,Z,varargin)
%KLEP_STEP2 Step 2.1 and 2.2 in [Klep et al.,2017]
%   V = KLEP_STEP2(n,m,Z) returns a tensor V of size m x n x m x n
%   (corresponding to a quadratic form f = sum v(i,k,j,l) xi xj yk yl) such
%   that f is a (2,2)-biform which is not divisible by any of the forms hi
%   lying in the kernel of the (random) matrix Z.
%
%   Main ideas: Derivative of 2x2 minors (Jacobian of Segre variety's
%   generators) and corresponding nullspace, symmetry nullspace, random
%   vector in intersection of both.
%
%   [V,K] = KLEP_STEP2(n,m,Z) also returns the involved (intersection of)
%   kernels.
%
%   V = KLEP_STEP2(n,m,Z,'example') simply returns V as explicitely given
%   in example 4.6 of [Klep et al., 2017].
%
%   See also KLEP_STEP1_1, KLEP_STEP1_2

d = n+m-2;
e = (n-1)*(m-1);

Z1 = Z(:,1:end-1);

N  = nchoosek(n,2);
M  = nchoosek(m,2);
ik = nchoosek(1:n,2); % list of (i,k) such that 1 <= i < k <= n, N possible
jl = nchoosek(1:m,2); % list of (j,l) such that 1 <= j < l <= m, M possible

% NB structure: vector [z_11, ..., z_1m, ..., z_nm] seen as matrix (z_ij)
% of size m x n (ith row = [z_1i, ..., z_ni] (rows and columns are
% switched)

% compute all possible combinations (i,k,j,l) with i<k and j<l
[A,B] = meshgrid(1:N,1:M);
ab    = [A(:),B(:)];
ikjl  = [ik(ab(:,1),:), jl(ab(:,2),:)]; % exahaustive list (i,k,j,l)

D_g = []; % Jacobian matrices
for t=1:N*M
	ik = ikjl(t,[1,2]); ki = flip(ik);
	jl = ikjl(t,[3,4]); lj = flip(jl);
	
	% {(Ipos,Jpos)} = {(i,j),(k,l),(i,l),(k,j)}, indices of variables
	% corresponding to non-zero entries in \nabla g_{ijkl} (stockage)
	Ipos = [ik,ik]; % indices of lines (position stockage)
	Jpos = [jl,lj]; % indices of columns (position stockage)
	
	% {(Ival,Jval)} = {(k,l),(i,j),(k,j),(i,l)}, indices of obtained
	% variables when differentiating wrt variables given by (Ival,Jval) in
	% expression z(i,j)z(k,l) - z(i,l)z(k,j)
	Ival = [ki,ki]; % indices of lines (position stocked)
	Jval = [lj,jl]; % indices of columns (position stocked)
	
	ivals   = sub2ind([m,n],Jval,Ival); % get corresponding values in Z = (z_ij) (NB: rows and lines switched)
	Dcoeffs = [1;1;-1;-1] .* Z1(ivals,:); % values of derivatives for g_t, at all the points z^(i), see Step 2.1
	
	% convert (m,n) matrix representation into (nm) vector
	ivec = sub2ind([m,n],Jpos,Ipos);
	
	% one column vector for each point 1,...,e
	ivec  = repmat(ivec,1,e);
	icols = repmat(1:e,4,1); % 4 because 4 nnz values (?)
	
	% gradient of one function gt for all points 1,...,e
	D_gt = sparse(ivec,icols,Dcoeffs(:),n*m,e);
	
	% sanity check against for loop
	%D_gt_2 = [];
	%for i=1:e %TODO: removing this loop should be possible, should it be done?
	%	D_gt_i = sparse(Ipos,Jpos,Dcoeffs(:,i),n,m); % sparse matrix representation for \nabla g_t(z^{(i)})
	%	D_gt_2 = cat(2,D_gt_2,D_gt_i(:));
	%	D_gt_i = []; % free memory, necessary?
	%end
	
	% final matrix: each column contains gradients of all gt, t=1,...,N, for
	% one point z^{(i)}.
	D_g = cat(1,D_g,D_gt);
end

% D_g matrix is of size (n*m*N*M) x e

% build kernel of gradient matrices [Grad g_1(Zi), ..., Grad g_M(Zi)]' and
% construct Kronecker matrices, see Step 2.2
ZW = zeros(n^2*m^2,d+1,e); % store the Kronecker products
for i=1:e
	% first build kernel basis
	Dg_i = reshape(D_g(:,i),n*m,N*M);
	K_i  = null(full(Dg_i)'); % basis (wi(1),...,wi(d+1)) see step 2.2 in paper by Klep
	
	% sanity check
	if size(K_i,2) ~= d+1
		error('kernel should be of dimension d+1 = %i, but is of dimension %i', d+1,size(K_i,2));
	end
	
	% compute kronecker product with data points, see step 2.3
	ZW(:,:,i) = kron(Z1(:,i),K_i); % omit z^{(e+1)}
end


%TODO: there must be a smarter way to implement this...
L = nchoosek(1:n*m,2); % list of indices 1 <= i < j <= n*m
%v = zeros(n*m,1);
I = sparse(1:n*m,1:n*m,1,n*m,n*m); % sparse representation of identity
E = [];
for i=1:size(L,1)
	l1 = L(i,1); l2 = L(i,2);
	E = cat(2,E,kron(I(:,l1),I(:,l2)) - kron(I(:,l2),I(:,l1)));
end

% choose vector in intersection of all kernels
ZWt     = permute(conj(ZW),[2,1,3]); % (hermitian) transpose each block
ZWt     = reshape(ZWt,e*(d+1),n^2*m^2); % up-down concatenation of transposed Kronecker matrices for all i, see step 2.2
Mt      = [ZWt;full(E)']; % concatenate with 
K_inter = null(Mt); % intersection of all kernels

% pick random vector in kernel
ga = rand(1,size(K_inter,2));
v  = sum(ga.*K_inter,2);
v  = reshape(v,[m,n,m,n]); %TODO: check if dimensions are correct!

if nargin > 3
	if n~=3 || m~=3
		error('example dimension mismatch: n,m = %i,%i (should be 3,3)',n,m);
	end
	% v comme dans l'exemple 4.6 du papier Klep
	v = [  5 ,  -2 ,  7/2,   2 ,  1/2,  7/4, -7/2, -3/4, -7/2; ...
		   -2 ,  -3 ,  -1 ,  1/2,  1/2, -5/2, -3/4,  5/2,  1/4; ...
		   7/2,  -1 ,   4 ,  7/4, -5/2,   1 , -7/2,  1/4,  -1 ; ...
		    2 ,  1/2,  7/4,  -2 , -1/2,  5/2,  -1 ,   0 , -5/4; ...
		   1/2,  1/2, -5/2, -1/2,   0 , -1/2,   0 ,  1/2,   0 ; ...
		   7/4, -5/2,   1 ,  5/2, -1/2,   1 , -5/4,   0 ,  -2 ; ...
		  -7/2, -3/4, -7/2,  -1 ,   0 , -5/4,   3 ,   1 ,   1 ; ...
		  -3/4,  5/2,  1/4,   0 ,  1/2,   0 ,   1 ,   2 ,   2 ; ...
		  -7/2,  1/4,  -1 , -5/4,   0 ,  -2 ,   1 ,   2 ,   2 ];
end

v = reshape(v,[m,n,m,n]); % TODO: check carefully correct dimensions!

end

