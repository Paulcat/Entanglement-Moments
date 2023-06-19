% A script to test the PnCP generation
% A state A is loaded
% A white Noise model is applied
% For growing value of noise, the entanglement of the state is tested in
% two ways : 1 - first, the program isentangled is called
%            2 - second, a PnCP is created, then applied to the state ;
%            check if the min eigenvalue is negative ; if it is, break

A = load('rho3x3.txt') ; 


	 % set options
	 options.mode     = 'hit-tol';
	 options.tol      = 1e-2;
	 options.ntest    = 30;
	 options.method   = 'klep';
	 options.maxorder = 2;
	 options.solver   = 'mosek';
	 options.toolbox  = 'yalmip';
	 options.verbose  = 0;



i=1;
while i>0.11
    i=i-0.1;
    B = WhiteNoise(A,i) ; 
    answer=isentangled([2,2],B,[3,3]);
    if answer==0
        ntemp = 0 ; 
        j = 10 ; % Max number of PnCP
        while ntemp > 0 || j > 0
            j = j -1;
            disp('PnCP generated to test entanglement')
            [phi,Deltas]=gen_PnCP(3,3,options);
            Mattemp=ApplyPnCP(B,phi);
            ntemp = min(eig(Mattemp));
        end
    else
        disp('PnCP generator not needed')
    end
end

if ntemp < 0
    disp('Entanglement detected thanks to PnCP generation')
else
    disp('Entanglement not detected through PnCP generation')
end

