function [lllll] =ketmat(n)
%created by Enky Oudot April 2018
% This programme will create a cube with dimension 2^N*2^N*N
% 

A=[0,0;1,1];
B=[0,1;0,1];
M=ones(2,2,2*n,n);


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


for k=1:n
lllll(:,:,k)=AA{k};
end

% on prend la difference de la somme



end



