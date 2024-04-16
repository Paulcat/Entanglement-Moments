function [PR,NPR] = PhysRegPoints1(N)


% Creates a grid in the Jx Jy Jz space of physical and unphysical points


%PAULI MATRICES
sigmax=[0 1; 1 0];
sigmay=[0 -1i; 1i 0];
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

        if indexs(ll, mm)==1               % indexs2 ??
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
%% 


%Preparation of the SDP

tot=10; %Number of points in each dimension for which we will run the SDP

%Values of Jx^2, Jy^2 and Jz^2 for which we sill run the SDP.
%Here we have to change the initial and final value depending on the number
%of qubits we have (if not maybe we won't be able to see the polytope...).
%Values here are for 3 qubits.
x0=linspace(0,2.25,tot);
y0=linspace(0,2.25,tot);
z0=linspace(0,2.25,tot);

%Feasibility problem.
object=1;
ops = sdpsettings('solver', 'mosek', 'verbose', 0);
rho = sdpvar(2^N,2^N) ;

PR=[]; %Points in the Physical Region
NPR=[]; %Points in the Not Physical Region


for ii=1:length(x0)
    
    for jj=1:length(y0)
        for kk=1:length(z0)
            
            %Constraints
            c=trace(rho)==1;
            c=c+(rho>=0);
            c=c+(trace(Jx*rho)==0); %Value of the first momenta aribtrarily fixed to 0.
            %c=c+(trace(Jy*rho)==0);
            c=c+(trace(Jz*rho)==0);
            c=c+(trace(J2x*rho)==x0(1,ii));
            c=c+(trace(J2y*rho)==y0(1,jj));
            c=c+(trace(J2z*rho)==z0(1,kk));

            %Run the SDP
            sol=optimize(c, object, ops);
            
            if sol.problem==0 %If the program finds a solution, the point is physical.
                
                PR=[PR; x0(1,ii), y0(1,jj), z0(1,kk)];
                
            else %If the program does not find a solution, the point is not physical.
                
                NPR=[NPR; x0(1,ii), y0(1,jj), z0(1,kk)];

                           
            end
        

            
            
            clear c
            
           
        end
        
    end
end


end

