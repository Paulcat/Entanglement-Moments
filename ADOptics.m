function M = ADOptics(rho, eta)

% Computes the effect of amplitude damping noise on a density matrix, according to a parameter eta
% The PartialTrace command needs the package QETLAB
% PROGRAM NOT WORKING

dims = size(rho, 1) ; % system dimension
dima = 2 ; % Ancilla dimension
[a, ad] = quantop(dims);
[b, bd] = quantop(dima);

Operator = expm(1i*pi*acos(eta)*(kron(ad,b)+kron(a,bd))) ; % Quantum operator of a beamsplitter that models effect of noise
Ancilla = zeros(dima);
Ancilla(1,1) = 1;


NoisyMatrix = transpose(conj(Operator))*kron(Ancilla, rho)*Operator ; % Transformation on the input state in the BS

M = PartialTrace(NoisyMatrix,1,2); % Resulting Noisy matrix ; I put Ancilla first, but it may be more natural to put it second and partial trace on the other subsystem, I just don't understand enough how PartialTrace function works. 
