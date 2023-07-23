function [Mat] = ApplyPnCP_left(dA,dB,rho,phi)
%APPLYPNCP_LEFT Apply (on the left) a PnCP on a given density matrix
%   Inputs:
%      - rho:   density matrix, rho := sum mijkl Eij x Ekl
%      - phi:   (PnCP) map
%      - dA,dB: bibpartite dimensions (rho must be of size dA*dB)
%
%   Output:
%      - Mat := sum mijkl phi(Eij) x Ekl


R = size(rho,1) ;  % assuming all matrices are square

if(dA*dB ~= R)
    warning("Dimensions mismatch between Ha, Hb, and Ha x Hb")
end

Mat = zeros(R);
for a=1:dA*dA
	Ea    = zeros(dA); Ea(a) = 1;
	
	for b=1:dB*dB
		Eb    = zeros(dB); Eb(b) = 1;
		Eab   = kron(Ea,Eb);

		% find coefficients of state wrt {Eab} basis
		val = trace(rho'*Eab) / norm(Eab,'fro')^2;
		
		% symmetrization
		Eas = (Ea+Ea')/2;

		% apply PnCP: acts on first (left) factor
		phiA = reshape(phi*Eas(:),[dA,dA]);
		Mat  = Mat + val * kron(phiA, Eb);
	end
end

end
