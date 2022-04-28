function  []= isentangled(order,state,dim)
% isentangled tries to check if a state is entangled using partial knowledge 
% of order order(1) for the first party and order(2) for the second party.

% state: a bipartite state with dimension dim(1) for the first party and dim(2) 
% dim: 1x2 vector which contain the dimension of the first and the second party
% order: 1x2 vector containing the order of the creation and anhilation operator one can measure for each party

c=moment(order,state,dim);

rho=sdpvar(size(state,1),size(state,2),'hermitian','complex');
%rho is a density matrix
 F2= [ rho >= 0] + [trace(rho) == 1];

%partial info of the moment of rho
F=[];
for i=1:size(c,1)
   F=F+ [trace(rho*c{i,1})==c{i,2}];
end
F3=PartialTranspose(rho,2,dim)>=0;
F=F+ F2+F3;

obj=1;
P = optproblem(F,obj,sdpsettings('solver','Mosek','verbose',0));
sol = minimize(P);
nbrmaps=0;
%this will be the size of the cellarraywiththemaps
g=value(rho);
if sol.problem == 0
 disp('we do not know if the state is entangled')
elseif sol.problem == 1
 disp('The state is entangled')
 for i=1:nbrmaps+1
           if i==1
           X(i)=min(eig(PartialTranspose(g,2,dim)));
              if X(1)<0
                  display('the state is NPT')
              end
           else
               X(i)=min(eig(somepncp(g,2,dim)));
               if X(i)<0
                  sprintf('Please enter the name for the %d xls file: ', i)
              end
           end
    
 end
 
 %add the test fot the other pncp maps
else
 disp('Something else happened')
end







end

