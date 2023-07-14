function [Mat] = ApplyPnCP_right(dA,dB,rho,phi)
%APPLYPNCP_RIGHT Apply (on the right) a PnCP on a given density matrix
%   Inputs:
%      - rho:   density matrix, rho := sum mijkl Eij x Ekl
%      - phi:   (PnCP) map
%      - dA,dB: bibpartite dimensions (rho must be of size dA*dB)
%
%   Output:
%      - Mat := sum mijkl Eij x phi(Ekl)


R = size(rho,1) ;  % assuming all matrices are square

if(dA*dB ~=R)
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
		
		% apply PnCP: acts on second factor
		phiB = reshape(phi*EBij(:),[dB,dB]);
		Mat = Mat + val * kron(EAij, phiB);
	end
end


% without Eij + Eji structure (wrong?)
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

