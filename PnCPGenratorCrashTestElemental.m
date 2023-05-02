% A script to test the PnCP generation to detect the entanglement of a
% state
% An entangled state A is loaded
% A PnCP is generated and applied ; the eigenvalues of the result are
% checked ; if at least one is negative, the entanglement is detected ;
% otherwise a new map is generated

rho = load('rho3x3.txt');
[dA,dB] = deal(3,3);
%rho = eye(9) ;
rho2 = WhiteNoise(rho,0.6);
phi_working = {};

nmax = 10;
times = zeros(40,1);
j = nmax ; % Max number of PnCP
while j > 0
    j = j-1;
    %disp('PnCP generated to test entanglement');

	 % set options
	 options.mode     = 'hit-tol';
	 options.tol      = 1e-2;
	 options.ntest    = 30;
	 options.method   = 'klep';
	 options.maxorder = 2;
	 options.solver   = 'mosek';
	 options.toolbox  = 'yalmip';
	 options.verbose  = 0;

	 tic;
    [phi,delta] = gen_PnCP(dA,dB,options);
	 tt = toc;
	 
	 % store computational time
	 times(nmax-j) = tt;
	 
    [~,Mattemp] = ApplyPnCP(dA,dB,rho2,phi);
    ntemp = min(eig(Mattemp));
    if ntemp < 0
        fprintf('%i. test: Entanglement detected\n',nmax-j);
		  phi_working{end+1} = phi;
          break
    else
        fprintf('%i. test: Not detected\n',nmax-j);
    end
end
