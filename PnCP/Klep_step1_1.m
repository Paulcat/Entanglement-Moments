function [Z,x,y] = Klep_step1_1(n,m,varargin)
%KLEP_STEP_1 Step 1.1 of algorithm described in [Klep et al.,2017]
%   [Z,X,Y] = KLEP_STEP1_1(n,m) returns random matrices
%
%       X = [x0, ..., xe],
%       Y = [y0, ..., ye], and
%       Z = [kron(x0,y0), .., kron(xe,ye)]
%
%   where e = (n-1)*(m-1) and each xi, yi are drawn randomly.
%
%   Z = KLEP_STEP_1(n,m,'example') simply returns X, Y and Z as explicitely
%   stated in example 4.6 of [Klep et al., 2017].
%
%   %   See also KLEP_STEP1_2, KLEP_STEP2

e = (n-1)*(m-1); % codimension

if nargin == 2 % random points
	x = rand(n,e+1);
	y = rand(m,e+1);
	
elseif nargin == 3 % example 4.6 in [Klep et al.]
	if n~=3 || m~=3
		error('example dimension mismatch: n,m = %i,%i (should be 3,3)',n,m);
	end
	x = [1 1 -1 1 2; 1 -1 1 1 -3; -1 1 1 1 3];
	y = [1 1 -1 1 -2; 1 -1 1 1 0; -1 1 1 1 2];
	
end

% form Kronecker products
Z = zeros(n*m,e+1);
for i=1:e+1
	Z(:,i) = kron(x(:,i),y(:,i));
end

