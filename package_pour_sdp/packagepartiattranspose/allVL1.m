function v = allVL1(n, L1, L1ops, MaxNbSol)
% All integer permutations with sum criteria
%
% function v=allVL1(n, L1); OR
% v=allVL1(n, L1, L1opt);
% v=allVL1(n, L1, L1opt, MaxNbSol);
% 
% INPUT
%    n: length of the vector
%    L1: target L1 norm
%    L1ops: optional string ('==' or '<=' or '<')
%           default value is '=='
%    MaxNbSol: integer, returns at most MaxNbSol permutations.
%    When MaxNbSol is NaN, allVL1 returns the total number of all possible
%    permutations, which is useful to check the feasibility before getting
%    the permutations.
% OUTPUT:
%    v: (m x n) array such as: sum(v,2) == L1,
%       (or <= or < depending on L1ops)                            
%       all elements of v is naturel numbers {0,1,...}
%       v contains all (=m) possible combinations
%       v is sorted by sum (L1 norm), then by dictionnary sorting criteria
%    class(v) is same as class(L1) 
% Algorithm:
%    Recursive
% Remark:
%    allVL1(n,L1-n)+1 for natural numbers defined as {1,2,...}
% Example:
%    This function can be used to generate all orders of all
%    multivariable polynomials of degree p in R^n:
%         Order = allVL1(n, p)
% Author: Bruno Luong
% Original, 30/nov/2007
% Version 1.1, 30/apr/2008: Add H1 line as suggested by John D'Errico
%         1.2, 17/may/2009: Possibility to get the number of permutations
%                           alone (set fourth parameter MaxNbSol to NaN)
%         1.3, 16/Sep/2009: Correct bug for number of solution
%         1.4, 18/Dec/2010: + non-recursive engine

global MaxCounter;

if nargin<3 || isempty(L1ops)
    L1ops = '==';
end

n = floor(n); % make sure n is integer

if n<1
    v = [];
    return
end

if nargin<4  || isempty(MaxNbSol)
    MaxCounter = Inf;
else
    MaxCounter = MaxNbSol;
end
Counter(0);

switch L1ops
    case {'==' '='},
        if isnan(MaxCounter)
            % return the number of solutions
            v = nchoosek(n+L1-1,L1); % nchoosek(n+L1-1,n-1)
        else
            v = allVL1eq(n, L1);
        end
    case '<=',% call allVL1eq for various sum targets
        if isnan(MaxCounter)
            % return the number of solutions
            %v = nchoosek(n+L1,L1)*factorial(n-L1); BUG <- 16/Sep/2009: 
            v = 0;
            for j=0:L1
                v = v + nchoosek(n+j-1,j);
            end
            % See pascal's 11th identity, the sum doesn't seem to
            % simplify to a fix formula
        else
            v = cell2mat(arrayfun(@(j) allVL1eq(n, j), (0:L1)', ...
                         'UniformOutput', false));
        end
    case '<',
        v = allVL1(n, L1-1, '<=', MaxCounter);
    otherwise
        error('allVL1: unknown L1ops')
end

end % allVL1