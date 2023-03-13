% A script to test the PnCP generation to detect the entanglement of a
% state
% An entangled state A is loaded
% A PnCP is generated and applied ; the eigenvalues of the result are
% checked ; if at least one is negative, the entanglement is detected ;
% otherwise a new map is detected

A = load('rho3x3.txt') ; 

j = 30 ; % Max number of PnCP
while j > 0
    j = j -1
    disp('PnCP generated to test entanglement')
    [phi,Deltas]=gen_PnCP(3,3,'hit-tol');
    Mattemp=ApplyPnCP(A,phi);
    ntemp = min(eig(Mattemp));
    if ntemp < 0
        disp('Entanglement detected thanks to PnCP generation')
    else
        disp('Not detected')
    end
end




