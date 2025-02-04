yalmip('clear');
X = sdpvar(2,2,'hermitian','complex');
Constraints = [X >= 0, trace(X) == 1];
Objective = -det(X);
options = sdpsettings('solver', 'mosek', 'verbose', 1);
diagnostics = optimize(Constraints, Objective, options);

if diagnostics.problem == 0
    disp('MOSEK works correctly for SDP.');
else
    disp(['MOSEK error: ', diagnostics.info]);
end
