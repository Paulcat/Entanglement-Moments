function [entvectors1,pptvectors1] = boxes5(h3)
%GET SOME SURFACE POINTS
[bf, surfacepoints]=boundaryFacets(h3);
%RUN SDP IN BOXES OF THE SURFACE POINTS
subtot=4; %New division. 

pptvectors1=[];
entvectors1=[];

for ii=1:length(surfacepoints(:,1))

    %Create boxes of new points around the surface points in order to gain
    %precision in the borders.
    %Epsilon was calculated in dmin (another script).
    x1=linspace(surfacepoints(ii,1)-epsilon, surfacepoints(ii,1)+epsilon,subtot);
    y1=linspace(surfacepoints(ii,2)-epsilon, surfacepoints(ii,2)+epsilon,subtot);
    z1=linspace(surfacepoints(ii,3)-epsilon, surfacepoints(ii,3)+epsilon,subtot);

    
    for ll=1:length(x1)
        for jj=1:length(y1)
            for kk=1:length(z1)
                if inShape(h2, x1(1,ll), y1(1,jj), z1(1,kk))==1 %avoid unphysical points
                    c=trace(rho)==1;
                    c=c+(rho>=0);
                    c=c+(trace(Jx*rho)==0); %value arbitrarily fixed to 0
                    % c=c+(trace(Jy*rho)==0);  
                    c=c+(trace(Jz*rho)==0);
                    c=c+(trace(J2x*rho)==x1(1,ll));
                    c=c+(trace(J2y*rho)==y1(1,jj));
                    c=c+(trace(J2z*rho)==z1(1,kk));
                    
                    %Partial transpose. Do all bipartitions (here it is
                    %just the case for 3 qubits).
                    c=c+(PartialTranspose(rho,1, [2,2,2])>=0); 
                    c=c+(PartialTranspose(rho,2, [2,2,2])>=0);
                    c=c+(PartialTranspose(rho,3, [2,2,2])>=0);
                    c=c+(PartialTranspose(rho,1, [2,4])>=0);
                    
                    %Run SDP
                    sol=optimize(c, object, ops);
                    if sol.problem==1
                        entvectors1=[entvectors1; x1(1,ll), y1(1,jj), z1(1,kk)];
                    else
                        pptvectors1=[pptvectors1; x1(1,ll), y1(1,jj), z1(1,kk)];
                    end
                    clear c
                else

                end
                
            end
        end
    
    end
    
    
end
end