% Test script

%% Test 1: generate PnCP from section 4.6 in [Klep et al.]

n = 3;
m = 3;
d = n+m-2;
%
N = nchoosek(n+1,2);
M = nchoosek(m+1,2);

% test random
[x,y,Z] = Klep_step1_1(n,m);
Vh      = Klep_step1_2(n,m,Z);
v       = Klep_step2(n,m,Z);

% % test example
% [x,y,Z] = Klep_step1_1(n,m,'example');
% Vh      = Klep_step1_2(n,m,Z,'example');
% v       = Klep_step2(n,m,Z,'example');

% check that we retrieve explicitely the polynomial obtained in [Klep...]
v     = reshape(v,n*m,n*m);
Vh    = reshape(Vh,n*m,d+1);
phi   = 2 * v + (Vh*Vh'); % coefficients of map
P_phi = zeros(M,N); % lines and columns switched (see comment in Klep_step_2)
E     = gen_sym_basis(n,m);
for i = 1:N
	for j = 1:M
		Eij        = E(:,:,i,j);
		P_phi(j,i) = trace(Eij'*phi);
	end
end

% solve sdp with Yalmip
[~,del,sol,Phi] = gen_one_PnCP(n,m,v,Vh);

% example of output
S = [0,0,1;0,0,0;1,0,0];
disp(Phi(m,S));


%% Test 2: generate a collection of PnCP (fixed initial random points)

n = 3;
m = 3;
d = n+m-2;
%
N = nchoosek(n+1,2);
M = nchoosek(m+1,2);

ntest = 10;

% specific initial points and resulting kernels
[~,~,Z]     = Klep_step1_1(n,m,'example');
[~,K,K1]    = Klep_step1_2(n,m,Z);
[~,K_inter] = Klep_step2(n,m,Z);

% dimensions
dK = size(K,2);
d1 = size(K1,2);
di = size(K_inter,2);
K  = reshape(K,[],1,dK);

% monitors
delta_list = zeros(ntest,1);
vf_list    = zeros(m,n,m,n,ntest);
vh_list    = zeros(m,n,d+1,ntest);

for a = 1:ntest
	fprintf('%i-',a);
	if ~mod(a,30)
		fprintf('\n');
	end
	% d random linear combination in K
	al = rand(1,d,dK);
	vj = sum(al.*K,3); % each column defines the linear form <vj',z>
	vj = reshape(vj,n,m,d); % row and columns are switched (cf Klep_step2.m)
	
	% one in K1
	be = rand(1,d1);
	v0 = sum(be.*K1,2);
	v0 = reshape(v0,m,n); % rows and columns are switched
	
	% all vectors
	vh = cat(3,v0,vj);
	
	% one in K_inter
	ga = rand(1,di);
	vf = sum(ga.*K_inter,2);
	vf = reshape(vf,[m,n,m,n]); % TODO: check if dimensions are in correct order!
	
	[Phi,delta] = gen_one_PnCP(n,m,vf,vh,'verbose',0);
	
	% store values
	delta_list(a)      = delta;
	vf_list(:,:,:,:,a) = vf;
	vh_list(:,:,:,a)   = vh;
end
fprintf('\n\n');

%% TODO: solve SDP with cvx instead of YALMIP

% going from map coefficients to polynomial coefficients
v = reshape(v,[n*m,n*m]);
E = gen_sym_basis(n,m);
C = zeros(M,N); % lines and columns switched (see comment in Klep_step_2)

for i = 1:N
	for j = 1:M
		Eij = E(:,:,i,j);
		
		C(j,i) = trace(Eij'*v);
	end
end

% and back
[~,En] = gen_sym_basis(n,m);

v2 = zeros(m*n,m*n);
for i=1:N
	for j=1:M
		Eij_n = En(:,:,i,j);
		v2   = v2 + C(j,i) * Eij_n;
	end
end
