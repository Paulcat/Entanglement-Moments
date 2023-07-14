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

[E,EA,EB] = gensymbasis(dA,dB);
Mat = zeros(R);
for i=1:size(E,1)
	for j=1:size(E,2)
		Eij = E{i,j};

		% corresponding tensorial factorization
		EAij = EA{i,j};
		EBij = EB{i,j};

		% % check
		%if norm(Eij-kron(EAij,EBij),'fro')/norm(Eij,'fro') > 1e-12
		%	error('incorrect matricial basis');
		%end

		% find coefficients of state in given basis
		val = trace(rho'*Eij);
		
		% apply PnCP: acts on first factor
		phiA = reshape(phi*EAij(:),[dA,dA]);
		Mat  = Mat + val * kron(phiA, EBij);
	end
end

end
