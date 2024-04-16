
%% States
%rho=load('rho3x3.txt');
%rho=WhiteNoise(delta,0.5);

%rho=rho2;
%rho=WhiteNoise(rho2,0.9);

%% Choice of subsystem


dA=3;
dB=3;


%% Looping through all the maps

s=zeros(20,1);
for i=1:20
    phi=PnCPs{i};

    PhiRho=ApplyPnCPSym_left(dA,dB,rho2,phi);
    s(i)=min(eig(PhiRho));
end
s