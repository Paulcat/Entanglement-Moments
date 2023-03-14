% A script to test the PnCP generation to detect the entanglement of a
% state
% An entangled state A is loaded
% A PnCP is generated and applied ; the eigenvalues of the result are
% checked ; if at least one is negative, the entanglement is detected ;
% otherwise a new map is detected

A = load('rho3x3.txt');

phi_working = {};

nmax = 100;
times = zeros(40,1);
j = nmax ; % Max number of PnCP
while j > 0
    j = j -1;
    %disp('PnCP generated to test entanglement');

	 % set options
	 options.mode = 'hit-tol';
	 options.tol = 1e-4;
	 options.ntest = 30;
	 options.method = 'klep';
	 options.maxorder = 2;
	 options.solver = 'mosek';
	 options.verbose = 0;

	 tic;
    [phi,delta] = gen_PnCP(3,3,options);
	 tt = toc;
	 
	 % store computational time
	 times(nmax-j) = tt;
	 
    Mattemp = ApplyPnCP(A,phi);
    ntemp = min(eig(Mattemp));
    if ntemp < 0
        fprintf('%i. test: Entanglement detected\n',nmax-j);
		  phi_working{end+1} = phi;
    else
        fprintf('%i. test: Not detected\n',nmax-j);
    end
end




