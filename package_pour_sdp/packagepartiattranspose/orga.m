function [L] = orga(N)
%   orga will select the elements of the Fock basis basiss which keep a coherence with the
%   qubit space plus the qubit space and put it in a matrix 
if isequal(N,2)
    L=[0,2;2,0];
else

for i=2:N

B{i-1}=allVL1(N, i, '==');
end
for k=1:N-1

kk=size(B{k});
s=1;
for i =1:kk(1)
l=0;
for j=1:kk(2)
if isequal(B{k}(i,j),1) || isequal(B{k}(i,j),0)
l=l+1;
end
end
if isequal(l,kk(2))

gff(s)=i;
s=s+1;
end
end

if isempty(gff)
gff=[];
end
pp=kk(1)*kk(2);
for i = 1:pp
index = true(1, size(B{k}, 1));
index(gff) = false;
y{k} = B{k}(index, :);
 
end
end
  L=[y{1};y{2}];
for l=1:size(y,2)-2
    L=[L;y{l+2}];
    
end
end
ff=ketmat(N);
for i=1:2^N
for j=1:N
lll(i,j)=ff(i,1,j);
end
end
L=[lll;L];
end

