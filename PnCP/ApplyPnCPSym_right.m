function [Mat] = ApplyPnCPSym_right(dA,dB,rho,phi)
%APPLYPNCPSYM_RIGHT Apply (extension of) a symmetric PnCP on the 2nd block
%   Inputs:
%      - rho:   density matrix, rho := sum mijkl Eij x Ekl
%      - phi:   (PnCP) map on symmetric matrices
%      - dA,dB: bipartite dimensions (rho must be of size dA*dB)
%
%   Output:
%      - Mat := sum mijkl Eij x phi((Ekl + Elk)/2)
%
%   Example:
%      If Phi is the transposition operator, and rho = rhoA x rhoB is a
%      (n,n) bipartite, separable state where rhoB is symmetric, then
%      APPLYPNCPSYM_RIGHT(n,n,rho,Phi) is the same as PartialTranspose(rho)
%      from QETLAB.
%
%   See also ApplyPnCPSym_left


R = size(rho,1) ;  % assuming all matrices are square

if(dA*dB ~=R)
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
		Ebs = (Eb+Eb')/2;
		
		% apply PnCP: acts on second (right) factor
		phiB = reshape(phi*Ebs(:),[dB,dB]);
		Mat = Mat + val * kron(Ea, phiB);
	end
end


% without Eij + Eji structure (wrong?--> right!)
% Mat = zeros(n);
% for i=1:dA
%     for j=1:dA
%         for k=1:dB
%             for l=1:dB
%                 EmatA=zeros(3); % HARD CODED
%                 EmatB=zeros(3); % HARD CODED
%                 EmatA(i,j)=1;
%                 EmatB(k,l)=1;
%                 Etemp = kron(EmatA,EmatB);
%                 scaltemp = trace(rho'*Etemp); % computing mijkl
%                 Mat = Mat + scaltemp*kron(EmatA,reshape(phi*EmatB(:),dB,dB)) ;% computes mijkl Eij * phi(Ekl)
%             end
%         end
%     end
% end

end
