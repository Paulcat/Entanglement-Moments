function [out2] = moment(order,rho,dim)
%moment(order,rho,dim) build a cell array out such that the first line of out
%contains all the operator a^i ad^j+h.c \otimes a^i1 ad^j1+h.c for
%j<=i<=order(1) and j1<=i1<=order(2). Dim is a vector such that dim(1)=dim Alice 
%and dim(2)=dim Bob
%The second line of out contains the expectation value Trace(rho*a^i ad^j+h.c) 



k=order(1);
k1=order(2);
N1=(dim(1));
N2=(dim(2));
if k>N1-1 || k1>N2-1
    display('order(1) and order(2) should be smaller than the local dimension -1')
    return
end
if N1*N2~=size(rho,1)
    display('The dimension of rho have to be equal to dim(1)*dim(2)')
    return
end
N1=N1-1;
N2=N2-1;
%we define the creation and anhilation operator on both side
a = diag(sqrt(1:N1),1); ad=a';
a2 = diag(sqrt(1:N2),1); ad2=a2';
%we define two arrays of indices  ind and ind5 for the local operator which grows with k
%and k1
ind1=repelem(0:k,1:k+1);
ind2=(0:((k+1)*(k+2)/2-1))-repelem((0:k).*((0:k)+1)/2,1:k+1);
ind=[ind1;ind2];
ind3=repelem(0:k1,1:k1+1);
ind4=(0:((k1+1)*(k1+2)/2-1))-repelem((0:k1).*((0:k1)+1)/2,1:k1+1);
ind5=[ind3;ind4];
% We build the cell array out where the first line contains the operator of
% Alice and the second, the operator of Bob
 out=cell(max([size(ind,2),size(ind5,2)]),2);
for i=1:max([size(ind,2),size(ind5,2)])
if i<=min(size(ind,2),size(ind5,2))
out{i,1}=((a^(ind(1,i)))*(ad^(ind(2,i))))+((a^(ind(1,i)))*(ad^(ind(2,i))))';
out{i,2}=((a2^(ind5(1,i)))*(ad2^(ind5(2,i))))+((a2^(ind5(1,i)))*(ad2^(ind5(2,i))))';
elseif size(ind,2)>=size(ind5,2)
out{i,1}=((a^(ind(1,i)))*(ad^(ind(2,i))))+((a^(ind(1,i)))*(ad^(ind(2,i))))';
else
out{i,2}=((a2^(ind5(1,i)))*(ad2^(ind5(2,i))))+((a2^(ind5(1,i)))*(ad2^(ind5(2,i))))';
end
end


k=1;
out2=cell(size(ind,2)*size(ind5,2),2);
      for i=1:(size(ind,2))
       for j=1:size(ind5,2)
          out2{k,1}=kron(out{i,1},out{j,2});
          out2{k,2}=trace(rho*out2{k,1});
    k=k+1;
    
       end
      end

end

