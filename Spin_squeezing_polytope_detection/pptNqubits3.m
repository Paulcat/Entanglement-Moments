function [h3,h4] = pptNqubits3(N,h2)



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

%% 



%Preparation of the SDP

tot=30; %Number of points in each dimension for which we will run the SDP

%Values of Jx^2, Jy^2 and Jz^2 for which we sill run the SDP.
%Here we have to change the initial and final value depending on the number
%of qubits we have (if not maybe we won't be able to see the polytope...).
%Values here are for 3 qubits.
x0=linspace(0.25,2.3,tot);
y0=linspace(0.25,2.3,tot);
z0=linspace(0.25,2.3,tot);

%Feasibility problem.
object=1;
ops = sdpsettings('solver', 'mosek', 'verbose', 0);
rho = sdpvar(2^N,2^N) ;

entvectors=[];
pptvectors=[];



for ii=1:length(x0)
    
    for jj=1:length(y0)
        for kk=1:length(z0)
            if inShape(h2, x0(1,ii), y0(1,jj), z0(1,kk))==1 %h2 is the physical region (calculated in another script). Here we just don't take into account the unphysical points.
                
                %Constraints.
                c=trace(rho)==1;
                c=c+(rho>=0);

                c=c+(trace(Jx*rho)==0);
                %disp(trace(Jy*rho)) ;
                %c=c+(trace(Jy*rho)==0);       % Because here trace(Jy*rho)
                %= 0 and matlab can not proceed constraints such as 0 ==0
                c=c+(trace(Jz*rho)==0);
                c=c+(trace(J2x*rho)==x0(1,ii));
                c=c+(trace(J2y*rho)==y0(1,jj));
                c=c+(trace(J2z*rho)==z0(1,kk));

                %PARTIAL TRANSPOSE. Here we should write the corresponding
                %partitions depending on the number of qubits (we have to
                %do the partial transpose with respect to every
                %partition!!! Here we have the example for 3 qubits.
%                 c=c+(PartialTranspose(rho,1, [2,2,2])>=0); 
%                 c=c+(PartialTranspose(rho,2, [2,2,2])>=0);
%                 c=c+(PartialTranspose(rho,3, [2,2,2])>=0);
%                 c=c+(PartialTranspose(rho,1, [2,4])>=0);
%                 c=c+(PartialTranspose(rho,1, [4,2])>=0);




                % For 4 qubits

                c=c+(PartialTranspose(rho,1, [2,2,2,2])>=0); 
                c=c+(PartialTranspose(rho,2, [2,2,2,2])>=0);
                c=c+(PartialTranspose(rho,3, [2,2,2,2])>=0);
                c=c+(PartialTranspose(rho,4, [2,2,2,2])>=0);

                c=c+(PartialTranspose(rho,1, [4,2,2])>=0);
                c=c+(PartialTranspose(rho,2, [4,2,2])>=0);
                c=c+(PartialTranspose(rho,3, [4,2,2])>=0);
                c=c+(PartialTranspose(rho,1, [2,4,2])>=0);
                c=c+(PartialTranspose(rho,2, [2,4,2])>=0);
                c=c+(PartialTranspose(rho,3, [2,4,2])>=0);
                c=c+(PartialTranspose(rho,1, [2,2,4])>=0);
                c=c+(PartialTranspose(rho,2, [2,2,4])>=0);
                c=c+(PartialTranspose(rho,3, [2,2,4])>=0);

                c=c+(PartialTranspose(rho,1, [2,8])>=0);
                c=c+(PartialTranspose(rho,1, [8,2])>=0);


                %Here we could add other PnCP maps, for example,
                %c=c+(BreuerMap(rho,[0,0])>=0); %The Breuer map.
                


                %Run the SDP
                sol=optimize(c, object, ops);
                            
                if sol.problem==0 %If the program finds a state, it is PPT
                    pptvectors=[pptvectors; x0(1,ii), y0(1,jj), z0(1,kk)]; 
               
                else %If the program is not able to find a state, it is entangled
                    entvectors=[entvectors; x0(1,ii), y0(1,jj), z0(1,kk)]; 
                    
                      
                end
            else
                                           
            end
                        
                                    
                                    
            clear c
            
           
         end
         
    end

h3=alphaShape(pptvectors, inf);
h4=alphaShape(entvectors, inf);


end

end