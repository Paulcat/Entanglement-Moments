function v = allVL1recurs(n, L1, head)
% function v=allVL1eq(n, L1);
% INPUT
%    n: length of the vector
%    L1: desired L1 norm
%    head: optional parameter to by concatenate in the first column
%          of the output
% OUTPUT:
%    if head is not defined
%      v: (m x n) array such as sum(v,2)==L1
%         all elements of v is naturel numbers {0,1,...}
%         v contains all (=m) possible combinations
%         v is (dictionnary) sorted
% Algorithm:
%    Recursive

global MaxCounter;

if n==1
    if Counter < MaxCounter
        v = L1;
    else
        v = zeros(0,1,class(L1));
    end
else % recursive call
    v = cell2mat(arrayfun(@(j) allVL1recurs(n-1, L1-j, j), (0:L1)', ...
                 'UniformOutput', false));
end

if nargin>=3 % add a head column
    v = [head+zeros(size(v,1),1,class(head)) v];
end

end % allVL1recurs