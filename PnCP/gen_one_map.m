function [phi,del,info,Phi] = gen_one_map(n,m,vf,vh,varargin)
%GEN_ONE_PNCP Generate a PnCP map
%   GEN_ONE_PNCP(n,m,vf,vh) solves a SDP program to generate a PnCP map.
%
%   Inputs:
%     -n : input dimension for map
%     -m : output dimension for map
%     -vf: quadratic form, seen as tensor, ie 
%             "vf(z x z)" = sum vf(i,k,j,l) z(i,j) z(k,l)
%           where z(i,j) = xi*yj (Segre variety)
%     -vh: sum-of-squares factors, seen as tensor, ie
%             "vh^2(zxz)" = sum (sum_a vh(i,j,a)*vh(k,l,a)) z(i,j) z(k,l)
%   
%   Options:
%     -'method'   : -'KLEP': solves the SDP hierarchy
%
%                        max  delta  s.t. (vh + delta*vf) * p0(l) = sos
%
%                    where p0 is a fixed sos polynomial of degree 2l (see
%                    step 3 of algorithm 4.3 in [Klep et al. 2017]).
%
%                   -'HILBERT': solves iteratively the feasibility programs
%
%                        min  1  s.t. (vh + delta0*vf) * sos1(l) = sos2
%
%                    where delta0 is decreased iteratively until success.
%     -'maxorder' : limit on the relaxation order (l)
%     -'tolerance': tolerance on delta (= must be greater than ...)
%     -'toolbox'  : 'YALMIP' | 'CVX' (not yet implemented)
%     -'solver'   : 'MOSEK' | 'SDPT3' | 'SeDuMi'
%     -'verbose'  :
%
%   Outputs:
%       - del: delta (or delta0) in the progrmas above
%       - phi: map corresponding to the form vh + del*vh. It is returned as
%              a matrix of size m*m x n*n (NB: there is a permutation
%              between v(i,k,j,l) (forms) and phi(i,j,k,l) (map))
%       - info:
%   
%   Needs YALMIP
%   Original author: Abhishek Bhardwaj (see functions step3*)
%
%   See also GEN_PNCP

d = n+m-2;

% default options
defaults = ...
	{'verbose', 2, ...
	 'maxorder', 2 ... % maximum relaxation order
	 'tolerance',1e-2,... % tolerance on delta
	 'method','klep', ...
	 'toolbox', 'yalmip', ...
	 'solver', 'mosek',...
	};

[verbose,maxorder,tol,method,toolbox,solver] = ...
	process_options(varargin,defaults{:});

switch toolbox
	case 'yalmip'
		if verbose
			fprintf('%5s Setting up YALMIP variables\n','');
		end
		% YALMIP solver
		x = sdpvar(n,1);
		y = sdpvar(m,1);
		z = kron(x,y);
		
		% quadratic forms
		vf = reshape(vf,n*m,n*m);
		vh = reshape(vh,n*m,d+1);
		f  = vf(:)' * kron(z,z); % positive not in <h0,...,hd> component
		h  = vh' * z; % h'*h SOS component

		opt = sdpsettings('solver',solver); % preferably mosek
		if verbose<=1
			opt.verbose = 0;
		else
			opt.verbose = 1;
		end

		switch method
			case 'klep'
				% coordinate norm relaxation see [Klep et al.]

				% sos perturbation
				sdpvar delta;

				% positive polynomial
				F = delta*f + (h'*h);

				% coordinate norm relaxation
				l = 1; % order of relaxation
				flag = 0;
				while l<=maxorder && ~flag
					if verbose
						fprintf('%5s Solving order %i of positive-non-sos relaxation problem\n','',l);
					end
					% Artin based relaxation with fixed denominator
					relax = F * (kron(x,y)'*kron(x,y))^l; % kron(x,y) is simply all monomials of bi-degree 1

					% Yalmip sos constraint
					constraint = sos(relax);

					% solve SDP
					[sol,u,Q,res] = solvesos(constraint,-delta,opt,delta);
					del = value(delta);

					% evaluate quality of solution
					flag = sol.problem==0 && floor(log10(res))<-5 && del > tol;
					
					l = l+1;
				end

			case 'hilbert'
				% general denominator + bisection method
				% method (3.1) in [Bhardwaj, 2020]

				% relax on degree of denominator
				l    = 1;
				flag = 0;
				while l <= maxorder && ~flag

					% monomials (for denominator)
					mon  = monolist([x;y],l);
					mon  = mon(2:end);
					nmon = size(mon,1);

					% symbolic Gram matrix for denominator
					G = sdpvar(nmon,nmon,'symmetric');
					D = mon'*G*mon;

					% initialize delta
					del = 1;

					% positive polynomial
					F = del*f + (h'*h);

					% sos constraint + trace constraint
					constraint = [sos(D*F); trace(G)==1];

					% bisection over del
					test = 0;
					while ~test && del>tol
						% solve SDP
						[sol,u,Q,res] = solvesos(constraint,[],opt,G(:));

						% quality of solution
						test = sol.problem==0 && floor(log10(res))<-6;
						del = del/2;
					end

					flag = test && del > tol;
					l = l+1;
				end

			case 'KKT'
				error('not implemented yet');
		end

		% info
		info.flag_sol = sol.problem;
		info.res      = res;
		info.success  = flag;

		% return PnCP map in correct format
		phi = del * vf  + (vh*vh');
		phi = reshape(phi,[m,n,m,n]); % check carefully correct dimensions
		phi = permute(phi,[1,3,2,4]);
		phi = reshape(phi,[m*m,n*n]);
		%
		Phi = @(m,S) reshape(phi*S(:),[m,m]);
		
		% verbose
		if verbose
			if sol.problem
				fprintf('%5s SDP was unsuccessfully solved: sol.problem: %i\n','',sol.problem);
			elseif ~flag
				fprintf('%5s SDP was succesfully solved, but map is not within tolerance: %d\n','',del);
			else
				fprintf('%5s Map within tolerance found: %d\n','',del);
			end
		end
				
		
	case 'cvx'
		error('not implemented yet');

		% sum-of-squares polynomial
		vsos = sum(reshape(vh,n*m,1,[]) .* reshape(vh,1,n*m,[]),3);
		vsos = reshape(vsos,m,n,m,n); % check carefully correct dimensions"
		
		l = 1; % order of relaxation
		
		% monomials matrices v(x,y)*v(x,y)'
		Hpos = genbihmat(n,m,1,1); % associated with the nonnegative polynomial
		Hsos = genbihmat(n,m,l,l); % associated with the sos denominator
		Hp   = genbihmat(n,m,l+1,l+1); % associated with the product of the two
		
		% sizes (for product)
		N = nchoosek(n+l,l+1);
		M = nchoosek(m+l,l+1);
		L = N*M;
		
		cvx_begin sdp
		variable Q(L,L) hermitian
		variable delta(1)
		%
		for i=1:size(Hp,1)
			for j=1:size(Hp,2)
				Hij = H{i,j}; % basis matrix
				convn(delta*vf + vsos, eye(3)) == trace(Hij'*Q);
			end
		end
end


end

