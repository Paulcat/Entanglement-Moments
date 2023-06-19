% A script to test the PnCP generation to detect the entanglement of a
% state
% An entangled state A is loaded
% A PnCP is generated and applied ; the eigenvalues of the result are
% checked ; if at least one is negative, the entanglement is detected ;
% otherwise a new map is generated

rho = load('rho3x3.txt');




[dA,dB] = deal(4,3);
%rho = eye(9) ;

% A = -1*ones(12,12);
% A(1,1)=5;
% A(1,2)=5;
% A(2,1)=5;
% A(2,2)=5;
% A(3,3)=5;
% A(3,12)=5;
% A(4,4)=5;
% A(4,10)=5;
% A(5,5)=5;
% A(5,6)=5;
% A(6,5)=5;
% A(6,6)=5;
% A(7,7)=5;
% A(7,9)=5;
% A(8,8)=5;
% A(8,11)=5;
% A(9,9)=5;
% A(9,7)=5;
% A(10,4)=5;
% A(10,10)=5;
% A(11,8)=5;
% A(11,11)=5;
% A(12,3)=5;
% A(12,12)=5;

rho2 = A/60;

rho2=rho;
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
	 times(nmax-j) = tt ;
	 
    [~,Mattemp] = ApplyPnCP(dA,dB,rho2,phi);
    ntemp = min(eig(Mattemp));
    if ntemp < 0
        fprintf('%i. test: Entanglement detected\n',nmax-j);
		  phi_working{end+1} = phi;
          break
    else
        fprintf('%i. test: Not detected\n',nmax-j);
    end
    phi
end
