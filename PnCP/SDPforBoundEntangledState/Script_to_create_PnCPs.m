	 % set options
	 options.mode     = 'hit-tol';
	 options.tol      = 1e-2;
	 options.ntest    = 30;
	 options.method   = 'klep';
	 options.maxorder = 2;
	 options.solver   = 'mosek';
	 options.toolbox  = 'yalmip';
	 options.verbose  = 0;


nbr_test=4;

Keeper_phi=zeros(16,16,4);

for i=1:nbr_test
    [phi,delta] = gen_PnCP(4,4,options);
    if phi~=0
        Keeper_phi(:,:,i)=phi;
    end
    i=i+1
end