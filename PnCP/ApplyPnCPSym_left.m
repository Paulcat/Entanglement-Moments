function [Mat] = ApplyPnCPSym_left(dA,dB,rho,phi)
%APPLYPNCPSYM_LEFT Apply (extension of) a symmetric PnCP on the 1st block
%   Inputs:
%      - rho:   density matrix, rho := sum mijkl Eij x Ekl
%      - phi:   (PnCP) map on symmetric matrices
%      - dA,dB: bipartite dimensions (rho must be of size dA*dB)
%
%   Output:
%      - Mat := sum mijkl phi((Eij+Eji)/2) x Ekl
%
%   Example:
%      If Phi is the transposition operator, and rho = rhoA x rhoB is a
%      (n,n) bipartite, separable state where rhoA is symmetric, then
%      APPLYPNCPSYM_LEFT(n,n,rho,Phi) is the same as PartialTranspose(rho)
%      from QETLAB.
%   
%   See also ApplyPnCPSym_right


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
