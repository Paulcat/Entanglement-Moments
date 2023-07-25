function [Mat] = ApplyPnCPSym_left(dA,dB,rho,phi)
%APPLYPNCPSYM_LEFT Apply (extension of) a symmetric PnCP on the 1st block.
%   M = APPLYPNCPSYM_LEFT(dA,dB,rho,phi) returns M = (PhiC x Id)(rho),
%   where PhiC is the extension to complex matrices of Phi:Sn -> Sm, which
%   is guaranteed to be PnCP if Phi is. No structure on rho is required.
%
%   Inputs:
%      - rho:   density matrix, rho := sum mijkl Eij x Ekl
%      - phi:   (PnCP) map on symmetric matrices
%      - dA,dB: bipartite dimensions (rho must be of size dA*dB)
%
%   Output:
%      - Mat := sum mijkl phi((Eij+Eji)/2) x Ekl
%
%   Warning:
%      Phi should not be a generic PnCP (eg the partial transpose): in that
%      case, APPLYPNCPSYM_LEFT actually applies ((Phi o Sym) x Id)(rho),
%      where Sym is the orthogonal projection on symmetric matrices. The
%      resulting map (Phi o Sym) might not be PnCP (?).
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
