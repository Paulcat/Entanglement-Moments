function test = test_positive(phi)
%TEST_POSITIVE Naive test for positivity of a map (no guarantee)
%   Test semidefinite positivity of Phi(A) for randomly generated sdp
%   matrices A


% number of test matrices
ntest = 100000;

% dimensions
[m2,n2] = size(phi); % m2 and n2 should be squares
m = sqrt(m2);
n = sqrt(n2);

i    = 0;
test = 1;


while i<ntest && test
	U = rand(n) + 1i*rand(n);
	A = U*U'; % sdp matrix
	
	B  = reshape(phi*A(:),[m,m]);
	Ev = eig(B);
	
	val = norm(imag(Ev),'fro') / norm(real(Ev),'fro');
	if val > 1e-15
		fprintf('image matrix may have complex eigenvalues: ||Re||/||Im|| = %d\n',val);
		test = 0;
	end
	
	if real(Ev(1)) < 0
		fprintf('image matrix may have negative eigenvalues: ev = %d\n',Ev(1));
		test = 0;
	end
	
	i = i+1;
end

end