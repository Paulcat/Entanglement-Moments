function [K] =matit(n)
%created by Enky Oudot April 2018
% This programme will create a matrix of one and 0 based on a phase averaging of the state (terms with 0 will cancel).
% The phase relation is define below and can be changed.
% The variable n is the number of modes.
% if d==1
% A=[0,0;1,1];
% B=[0,1;0,1];
% M=ones(2,2,2*n,n);
% elseif d==2
% A=[0,0,0;1,1];
% B=[0,1;0,1];


A=[0,0;1,1];
B=[0,1;0,1];
M=ones(2,2,2*n,n);
if isequal(n,2)
  K=[1,0,0,0;
      0,1,1,0;
      0,1,1,0;
      0,0,0,1];  
else

for i=1:n
for j=1:n
    if j==i 
M(:,:,i,j)=A;
    end
end
end
for i=n+1:2*n
for j=1:n
    if (j+n)==i 
M(:,:,i,j)=B;
    end
end
end

for i=1:2*n
    AA{i}=kron(M(:,:,i,1),M(:,:,i,2));
    for j=1:n-2
AA{i}=kron(AA{i},M(:,:,i,j+2));
    end
end
% on a donc les lignes dans les n premieres matrices et les colonnes dans
% les n suivantes on choisit ici de sommer les lignes et les collones (les ket et les bras)
for k=1:n
lllll(:,:,k)=AA{k};
end
for k=n+1:2*n
lllllll(:,:,k-n)=AA{k};
end

% on prend la difference de la somme(selon la relation de phase (ici on veut la somme))

%S=sum(lllll,3)-sum(lllllll,3);

% ket=sum(lllll,3)+lllll(:,:,2)+lllll(:,:,3)+lllll(:,:,4);
% bra=lllllll(:,:,1)+lllllll(:,:,2)+lllllll(:,:,3)+lllll(:,:,4);
ket=sum(lllll,3);
bra=sum(lllllll,3);
S=ket-bra;
for i=1:size(S,1)
for j=1:size(S,2)
   if isequal(S(i,j),0)
       K(i,j)=1;
   else
       K(i,j)=0;
   end
end
end

end
end
