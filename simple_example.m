%basic example with a Bellstate
psi=[1/sqrt(2),0,0,1/sqrt(2)];
state=psi'*psi;
isentangled([1,],state,[2,2])
%add other example