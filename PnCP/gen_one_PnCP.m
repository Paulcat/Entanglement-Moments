function [phi,del,sol,Phi] = gen_one_PnCP(n,m,vf,vh,varargin)
%GEN_ONE_PNCP Generate a PnCP map
%   Needs YALMIP
%   Original author: Abhishek Bhardwaj
%
%   Inputs:
%       - n, m: dimensions
%       - vf:   quadratic form, seen as tensor, ie 
%               "vf(z x z)" = sum vf(i,j,k,l) z(i,j) z(k,l)
%               where z(i,j) = x(i)*y(j) (Segre variety)
%       - vh:   sum-of-squares factors, seen as tensor, ie
%               "vh^2(zxz)" = sum (sum_a vh(i,j,a)*vh(k,l,a)) z(i,j) z(k,l)
%
%   Options:
%       - verbose
%       - toolbox
%       - solver

d = n+m-2;

% default options
defaults = ...
	{'verbose', 2, ...
	 'order-max', 1 ... % maximum relaxation order (inutile)
	 'toolbox', 'yalmip', ...
	 'solver', 'mosek'
	};

[verbose,order,toolbox,solver] = process_options(varargin,defaults{:});



switch toolbox
	case 'yalmip'
		% YALMIP solver
		x = sdpvar(n,1);
		y = sdpvar(m,1);
		z = kron(x,y);
		
		
		% quadratic forms
		vf = reshape(vf,n*m,n*m);
		vh = reshape(vh,n*m,d+1);
		f  = vf(:)' * kron(z,z); % positive not in <h0,...,hd> component
		h  = vh' * z; % h'*h SOS component
		
		sdpvar delta
		
		mon = monolist([x;y],d);
		mon = mon(2:end);
		
		% positive polynomial
		F = delta * f + (h'*h);
		
		% Artin-based relaxation
		l = 1; % order of relaxation
		relax = F * (kron(x,y)'*kron(x,y))^l;
		
		% Yalmip sos constraint
		constraint = sos(relax);
		
		% sdp solver
		opt = sdpsettings('solver',solver); % preferably mosek
		opt.verbose = verbose;
		%
		[sol,u,Q,res] = solvesos(constraint,-delta,opt,delta);
		del = value(delta);
		
		% TODO: implement a loop if tolerance not reached to go to
		% higher relaxations orders
		
		% return PnCP map in correct format
		phi = del * vf  + (vh*vh');
		%phi = 2 * vf + (vh*vh');
		phi = reshape(phi,[m,n,m,n]); % check carefully correct dimensions"
		phi = permute(phi,[1,3,2,4]);
		phi = reshape(phi,[m*m,n*n]);
		%
		Phi = @(m,S) reshape(phi*S(:),[m,m]);
		
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


function varargout = process_options(params,varargin)
% A simplified version of the parseparams code in octave
% author:
% Jean-Francois Cardoso
%
% original authors:
% Alexander Barth
% Aida Alvera

% default options
names  = varargin(1:2:end);
defvalues = varargin(2:2:end);

% specified options
pnames  = params(1:2:end);
values = params(2:2:end);
if (length(pnames)~=length(values)) || ~iscellstr(pnames)
	error('options must be given as name-value pairs');
end

% match
varargout = defvalues;
for i = 1:length(pnames)
	pname = pnames{i};
	id = find(strcmp(names,pname));

	% set values
	if id==0
		error('unknown option: %s', pname);
	end
	varargout{id} = values{i};
end

end

