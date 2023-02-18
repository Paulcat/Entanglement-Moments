function [Phi,del,sol] = gen_PnCP(n,m,vf,vh,varargin)
%GEN_PNCP Generate a PnCP map
%   Needs YALMIP
%   Copyright: Abhishek Bhardwaj

d = n+m-2;



% process input options % TODO: put in separate function
defaults = {'verbose',2};
% default options
names  = defaults(1:2:end);
values = defaults(2:2:end);
% specified options
snames  = varargin(1:2:end);
svalues = varargin(2:2:end);

for a = 1:length(snames)
	sname = snames{a};
	id = find(strcmp(names,sname));
	
	% set values
	if id==0
		error('unknown option: %s', sname);
	end
	values{id} = svalues{id};
end
verb = values{1}; % hack, for now


% YALMIP solver
x = sdpvar(n,1);
y = sdpvar(m,1);
z = kron(x,y);

sdpvar delta

mon = monolist([x;y],d);
mon = mon(2:end);

% quadratic forms
vf = reshape(vf,n*m,n*m);
vh = reshape(vh,n*m,d+1);
f  = vf(:)' * kron(z,z); % positive not in <h0,...,hd> component
h  = vh' * z; % h'*h SOS component

F = delta * f + (h'*h);

% Artin-based relaxation
l = 1; % order of relaxation
relax = F * (kron(x,y)'*kron(x,y))^l;

% Yalmip sos constraint
constraint = sos(relax);

% sdp solver
opt = sdpsettings('solver','sdpt3'); % preferably mosek
opt.verbose = verb;
%
[sol,u,Q,res] = solvesos(constraint,-delta,opt,delta);
del = value(delta);

% return PnCP map in correct format
phi = del * vf  + (vh*vh');
%phi = 2 * vf + (vh*vh');
phi = reshape(phi,[m,n,m,n]); % check carefully correct dimensions"
phi = permute(phi,[1,3,2,4]);
phi = reshape(phi,[m*m,n*n]);
%
Phi = @(m,S) reshape(phi*S(:),[m,m]);


end

