function [V,K,K1] = Klep_step1_2(n,m,Z,varargin)
%KLEP_STEP1_2  Step 1.2 and 1.3 of algorithm described in [Klep et al.]
%		Find vectors in kernel, define linear forms:
%     Returns a tensor V of size m x n x (m+n-2) containing all the vectors
%     vi of size m*n (seen as matrices m x n) representing the linear forms
%     hi

d = n+m-2;
e = (n-1)*(m-1);
if size(Z,2) ~= e+1
	error('dimension mismatch: %i columns in Z, should be %i', size(Z,2),e+1);
end

% Step 1.2
% kernel associated with our random points
K  = null(Z'); % could be implemented with svd
dK = size(K,2); %dimension of kernel. TODO: apparently dK = d, is that always true?

% take random linear combination of vectors in K
K2 = reshape(K,n*m,1,dK);
al = rand(1,d,dK); % d random coefficients
al = al./sum(al,3); % normalize
Vj = sum(al.*K2,3); % each resulting column is a random convex combination of all the columns in K
Vj = reshape(Vj,m,n,d); % row and columns are switched in my encoding

if nargin > 3 % example 4.6 in [Klep et al.]
	if n~=3 || m~=3
		error('example dimension mismatch: n,m = %i,%i (should be 3,3)',n,m);
	end
	Vj(:,:,1) = reshape([0;2;3;-2;3;0;-3;0;-3],m,n);
	Vj(:,:,2) = reshape([-3;7;0;-7;-3;1;0;-1;6],m,n);
	Vj(:,:,3) = reshape([9;-14;0;14;-3;2;0;-2;-6],m,n);
	Vj(:,:,4) = reshape([0;6;0;-6;-6;0;0;0;6],m,n);
end

%warning('We must make sure that Ker(Vj*) intersects sufficiently with Segre variety, see Step 1.2');

% Step 1.3
% kernel associated with our random points minus last
Z1 = Z(:,1:end-1);
K1 = null(Z1');
d1 = size(K1,2);

% take one random vector in K1
be = rand(1,d1);
v0 = sum(be.*K1,2);
v0 = reshape(v0,m,n);

if nargin > 3
	v0 = [-2,-2,1;2,0,0;-1,0,2];
end

%warning('same test to perform, see Step 1.3');

V = cat(3,v0,Vj); % V of size m x n x (d+1), contains all the coefficietns of h0,...,hd

