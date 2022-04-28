function rho = RandomDensityMatrix(dim,varargin)
 
% set optional argument defaults: re=0, k=dim, dist='hs'
[re,k,dist] = opt_args({ 0, dim, 'haar' },varargin{:});

% Haar/Hilbert-Schmidt measure
gin = randn(dim,k);
if(~re)
    gin = gin + 1i*randn(dim,k);
end
if(strcmpi(dist,'bures')) % Bures measure
    gin = (RandomUnitary(dim,re) + eye(dim))*gin;
end
 
rho = gin*gin';
rho = rho/trace(rho);
