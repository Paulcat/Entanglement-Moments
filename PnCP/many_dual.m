function [a,b,c,d]=many_dual(N, x0, y0, z0, maps)

% (x0, y0, z0) is the point in the (Jx^2, Jy^2, Jz^2) space that we want to detect.

%(a, b, c) is the normal vector of the plane.

%d is the constant of the plane so that a*x+b*y+c*z=d defines a plane (our entanglement witness).

%PAULI MATRICES
sigmax=[0 1; 1 0];
sigmay=[0 -i; i 0];
sigmaz=[1 0; 0 -1];


%FIRST ORDER MOMENTA

indexs=[1]; %This is used to make the first order momenta

for ss=1:(N-1)
    indexs=[indexs,0];
end

indexs=unique(perms(indexs), 'rows'); %do all permutations possible


Delta=zeros(2,2,N);
Zeta=zeros(2,2,N);
Eta=zeros(2,2,N);


Jx=zeros(2^N);
Jy=zeros(2^N);
Jz=zeros(2^N);


for ll=1:N
    
    for mm=1:N

        if indexs(ll, mm)==1   %indexs2
            Delta(:,:,mm)=sigmax;
            Zeta(:,:,mm)=sigmay;
            Eta(:,:,mm)=sigmaz;

            
        else
            Delta(:,:,mm)=eye(2);
            Zeta(:,:,mm)=eye(2);
            Eta(:,:,mm)=eye(2);
        end
    end
    
    Rx=Delta(:, :, 1);
    Ry=Zeta(:, :, 1);
    Rz=Eta(:, :, 1);
    for nn=2:N
        Rx=kron(Rx, Delta(:, :, nn));
        Ry=kron(Ry, Zeta(:, :, nn));
        Rz=kron(Rz, Eta(:, :, nn));
        
    end

    Jx=Jx+1/2.*Rx;
    Jy=Jy+1/2.*Ry;
    Jz=Jz+1/2.*Rz;
end

%SECOND ORDER MOMENTA

J2x=Jx*Jx;
J2y=Jy*Jy;
J2z=Jz*Jz;





object=1;
ops = sdpsettings('solver', 'mosek', 'verbose', 0);
rho = sdpvar(2^N,2^N) ;

%Constraints.

%Here we add other PnCP maps, for example,
F=[trace(J2x*rho)==x0, trace(J2y*rho)==y0, trace(J2z*rho)==z0, trace(rho)==1, rho>=0, trace(Jx*rho)==0, trace(Jz*rho)==0];
for pp=1:length(maps)
    F=[F,(ApplyPnCP_right(2,4,rho,maps{1, pp})>=0)];
    F=[F, (ApplyPnCP_left(4,2,rho,maps{1, pp})>=0)];
end

sol=optimize(F, object, ops);

a=dual(F(1));
b=dual(F(2));
c=dual(F(3));


rhod = sdpvar(2^N,2^N) ;
objd=a*trace(J2x*rhod)+b*trace(J2y*rhod)+c*trace(J2z*rhod);

Fd=[trace(rhod)==1, rhod>=0, trace(Jx*rhod)==0, trace(Jz*rhod)==0];

for jj=1:length(maps)
    Fd=[Fd,(ApplyPnCP_right(2,4,rhod,maps{1, jj})>=0)];
    Fd=[Fd, (ApplyPnCP_left(4,2,rhod,maps{1, jj})>=0)];
end


%[cd, objd, primals]=dualize(c, object);

sold=optimize(Fd, objd, ops);

dm=value(rhod);
d=a*trace(J2x*dm)+b*trace(J2y*dm)+c*trace(J2z*dm);

end