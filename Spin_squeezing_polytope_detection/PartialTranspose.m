function Xpt = PartialTranspose(X,varargin)

dX = size(X);
sdX = round(sqrt(dX));

% set optional argument defaults: sys=2, dim=round(sqrt(length(X)))
[sys,dim] = opt_args({ 2, [sdX(1) sdX(1);sdX(2) sdX(2)] },varargin{:});

num_sys = length(dim);

% allow the user to enter a single number for dim
if(num_sys == 1)
    dim = [dim,dX(1)/dim];
    if abs(dim(2) - round(dim(2))) >= 2*dX(1)*eps
        error('PartialTranspose:InvalidDim','If DIM is a scalar, X must be square and DIM must evenly divide length(X); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
    num_sys = 2;
end

% allow the user to enter a vector for dim if X is square
if(min(size(dim)) == 1)
    dim = dim(:).'; % force dim to be a row vector
    dim = [dim;dim];
end

% prepare the partial transposition
prod_dimR = prod(dim(1,:));
prod_dimC = prod(dim(2,:));
sub_prodR = prod(dim(1,sys));
sub_prodC = prod(dim(2,sys));
sub_sys_vecR = prod_dimR*ones(1,sub_prodR)/sub_prodR;
sub_sys_vecC = prod_dimC*ones(1,sub_prodC)/sub_prodC;
perm = [sys,setdiff(1:num_sys,sys)];

Xpt = PermuteSystems(X,perm,dim); % permute the subsystems so that we just have to do the partial transpose on the first (potentially larger) subsystem

if(isnumeric(Xpt)) % if the input is a numeric matrix, perform the partial transpose operation the fastest way we know how
    Xpt = cell2mat(mat2cell(Xpt, sub_sys_vecR, sub_sys_vecC).'); % partial transpose on first subsystem
else % if the input is not numeric (such as a variable in a semidefinite program), do a slower method that avoids mat2cell (mat2cell doesn't like non-numeric arrays)
    Xpt = reshape(permute(reshape(Xpt,[sub_sys_vecR(1),sub_prodR,sub_sys_vecC(1),sub_prodC]),[1,4,3,2]),[sub_sys_vecR(1)*sub_prodC,sub_sys_vecC(1)*sub_prodR]);
end

% return the subsystems back to their original positions
dim(:,sys) = dim([2,1],sys);
dim = dim(:,perm);
Xpt = PermuteSystems(Xpt,perm,dim,0,1);
