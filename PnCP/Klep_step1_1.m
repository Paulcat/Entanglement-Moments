function [x,y,Z] = Klep_step1_1(n,m,varargin)
%KLEP_STEP_1 Step 1.1 of algorithm described in [Klep et al.]

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

