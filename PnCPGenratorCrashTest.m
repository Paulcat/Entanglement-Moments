% A script to test the PnCP generation
% A state A is loaded
% A white Noise moel is applied
% For growing value of noise, the entanglement of the state is tested in
% two ways : 1 - first, the program isentangled is called
%            2 - second, a PnCP is created, then applied to the state ;
%            check if the min eigenvalue is negative ; if it is, break

A = load('rho3x3.txt') ; 

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
            [phi,Deltas]=gen_PnCP(3,3,'hit-tol');
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

