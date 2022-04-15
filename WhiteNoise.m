function M = WhiteNoise(rho, eta)
%Computes the effect of white noise on a density matrix, according to a parameter eta

dim = size(rho, 1) ; % size of input density matrix
M = eta*rho + ((1-eta)/dim)*eye(dim) ; % Noise channel