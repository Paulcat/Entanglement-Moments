function [h2] = PhysRegPointsBoxes2(N,PR,NPR,epsilon)

subtot=4;

PR1=[];
NPR1=[];

rho = sdpvar(2^N,2^N) ;

% Enhances the precision of PhysRegPoints1

for ii=1:length(PR(:,1))
    for mm=1:length(NPR(:,1))
        dis=distance(PR(ii,1), PR(ii,2), PR(ii,3), NPR(mm,1), NPR(mm,2), NPR(mm,3));
        if dis<epsilon
            
            %Do some boxes around the physical points closer to not physical states
            %in order to gain precision
            x1=linspace(PR(ii,1)-epsilon, PR(ii,1)+epsilon,subtot);
            y1=linspace(PR(ii,2)-epsilon, PR(ii,2)+epsilon,subtot);
            z1=linspace(PR(ii,3)-epsilon, PR(ii,3)+epsilon,subtot);

            for ll=1:length(x1)
                for jj=1:length(y1)
                    for kk=1:length(z1)
                        c=trace(rho)==1;
                        c=c+(rho>=0);
                        
                        c=c+(trace(Jx*rho)==0);
                        c=c+(trace(Jy*rho)==0);
                        c=c+(trace(Jz*rho)==0);
                        c=c+(trace(J2x*rho)==x1(1,ll));
                        c=c+(trace(J2y*rho)==y1(1,jj));
                        c=c+(trace(J2z*rho)==z1(1,kk));
                        
                        sol=optimize(c, object, ops);
                        if sol.problem==1
                            
                            NPR1=[NPR1; x1(1,ll), y1(1,jj), z1(1,kk)];
                        else
                            
                            PR1=[PR1; x1(1,ll), y1(1,jj), z1(1,kk)];
                        end
                        clear c
                    end
                end
            
            end
        else
            a=0;
        
        end
    end
end


PRtot=[PR; PR1]; %Concatenate all the points

h2=alphaShape(PRtot, inf); %Create the volume (it is a cone) of physical states

end
