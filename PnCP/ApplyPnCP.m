function [Mat] = ApplyPnCP(rho,phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply a PnCP on a given density matrix
%   rho = mijkl Eij x Ekl
%   Mat = mijkl Eij x phi(Ekl)
%   rho := density matrix
%   phi := PnCP map
n = size(rho,1) ;  % assuming all matrices are square

% TO DO : AUTOMATIZE SIZE OF da and db

da = 3;
db = 3;
if(da*db ~=n)
    warning("Dimensions mismatch between Ha, Hb, and Ha x Hb")
end

% TO DO : REDUCE BASIS OF EIJ
% TO DO : use multi indices

Mat = zeros(n);
for i=1:da
    for j=1:da
        for k=1:db
            for l=1:db
                EmatA=zeros(3); % HARD CODED
                EmatB=zeros(3); % HARD CODED
                EmatA(i,j)=1;
                EmatB(k,l)=1;
                Etemp = kron(EmatA,EmatB) ; 
                scaltemp = trace(rho'*Etemp); % computing mijkl
                Mat = Mat + scaltemp*kron(EmatA,reshape(phi*EmatB(:),db,db)) ;% computes mijkl Eij * phi(Ekl)
            end
        end
    end
end

end

