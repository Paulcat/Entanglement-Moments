rho=sdpvar(8,8,'hermitian', 'complex');

sigmax=[0 1; 1 0];
sigmay=[0 -i; i 0];
sigmaz=[1 0; 0 -1];

N=3; %3 qubits:



indexs2=[0,0,1];

indexs2=unique(perms(indexs2), 'rows');


Delta=zeros(2,2,N);
Zeta=zeros(2,2,N);
Eta=zeros(2,2,N);



Jx=zeros(8);
Jy=zeros(8);
Jz=zeros(8);




%Do the same for Jx Jy Jz (without the square)


for ll=1:N
    
    for mm=1:N

        if indexs2(ll, mm)==1
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

J2x=Jx*Jx;
J2y=Jy*Jy;
J2z=Jz*Jz;


%Initial conditions:
tot=17;


y0=linspace(0.25,2.3,tot);
z0=linspace(0.25,2.3,tot);
x0=linspace(0.25,2.3,tot);

object=1;
ops = sdpsettings('solver', 'mosek', 'verbose', 0);

breuer3=[];
notbreuer3=[];
numericalerrorpoint=[];
EIGEN=[];

for ii=1:length(x0)
    
    for jj=1:length(y0)
        for kk=1:length(z0)
            if inShape(h2, x0(1,ii), y0(1,jj), z0(1,kk))==1 %h2 is the physical region (calculated in another script)
                c=trace(rho)==1;
                c=c+(rho>=0);
                c=c+(BreuerMap(rho,R,[0,0])>=0); 
                
                c=c+(trace(Jx*rho)==0);
                c=c+(trace(Jy*rho)==0);
                c=c+(trace(Jz*rho)==0);
                c=c+(trace(J2x*rho)== x0(1,ii));
                c=c+(trace(J2y*rho)== y0(1,jj));
                c=c+(trace(J2z*rho)== z0(1,kk));
                sol=optimize(c, object, ops);

                
                if sol.problem==0
                    dm=value(rho);
                    eigenvalues=eigs(dm);
                    if max(abs(eigenvalues))>10^(-6)

                        breuer3=[breuer3; x0(1,ii), y0(1,jj), z0(1,kk)]; 
                        EIGEN=[EIGEN; eigenvalues'];

                    else
                        numericalerrorpoint=[numericalerrorpoint; x0(1,ii), y0(1,jj), z0(1,kk)];
                    end

                    


                else
                   notbreuer3=[notbreuer3; x0(1,ii), y0(1,jj), z0(1,kk)]; %PRtot(ii,1), PRtot(ii,2), PRtot(ii,3)]; %
                    
                      
                end
            else
                a=0;                             
            end
                        
                                    
                                    
            clear c
            
           
         end
         
     end
end

scatter3(breuer3(:,1), breuer3(:,2), breuer3(:,3), 70, '.', 'MarkerEdgeColor', 'yellow')
hold on
scatter3(notbreuer3(:,1), notbreuer3(:,2), notbreuer3(:,3), 70, '.', 'MarkerEdgeColor', 'cyan')