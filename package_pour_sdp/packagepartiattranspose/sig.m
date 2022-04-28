function [BBBB] = sig(n,k,l,alpha1,alpha2)
%UNTITLED Summary of this function goes here
%   this function gives the observable sigma_alpha1 otimes sigma_alpha2 for the mode k and
%   l in a n qubit space.
A=dec2bin(2^n-1:-1:0)-'0';
A=flip(A);

 for i=1:size(A,1)
for j=1:size(A,1)
BBB(i,j,:)=[A(i,:),A(j,:)];
end
 end



t=1;
for i=1:n
if not(i==k) && not(i==l)
ff(t)=i;
t=t+1;
end
end
for i=1:size(A,1)
for j=1:size(A,1)
if isequal(BBB(i,j,ff(1)),BBB(i,j,ff(1)+n)) &&  isequal(BBB(i,j,ff(2)),BBB(i,j,ff(2)+n))
    if isequal(BBB(i,j,k),BBB(i,j,k+n)) && isequal(BBB(i,j,l),BBB(i,j,l+n))
c1=1;
    else
        c1=0;
    end
BBBB(i,j)=4*exp(-alpha1^2)*exp(-alpha2^2)*alpha1^(BBB(i,j,k))*alpha1^(BBB(i,j,k+4))*alpha2^(BBB(i,j,l))*alpha2^(BBB(i,j,l+4))+c1;
end
end
end
 

for i=1:size(A,1)
for j=1:size(A,1)
if isequal(BBB(i,j,ff(1)),BBB(i,j,ff(1)+4)) &&  isequal(BBB(i,j,ff(2)),BBB(i,j,ff(2)+4))
     if  isequal(BBB(i,j,l),BBB(i,j,l+4))

BBBB(i,j)=BBBB(i,j)-2*exp(-alpha1^2)*alpha1^(BBB(i,j,k))*alpha1^(BBB(i,j,k+4));

     end


end
end
end



for i=1:size(A,1)
for j=1:size(A,1)
if isequal(BBB(i,j,ff(1)),BBB(i,j,ff(1)+4)) &&  isequal(BBB(i,j,ff(2)),BBB(i,j,ff(2)+4))
     if  isequal(BBB(i,j,k),BBB(i,j,k+4))

BBBB(i,j)=BBBB(i,j)-2*exp(-alpha2^2)*alpha2^(BBB(i,j,l))*alpha2^(BBB(i,j,l+4));

     end


end
end
end


 
 
 
end

