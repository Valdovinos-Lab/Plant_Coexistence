% Code created for simulating invasion of pollinator species into
% plant-pollinator networks, using Valdovinos et al. 2013, 2016 model
% Authors: Fernanda S. Valdovinos & Pablo Moisset de Espanes
% Last Modification: May 5, 2018
% Published: Nature Communications 2018 (see paper for more in detail
% explanation of parameters, variables, and simulation design)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Valdovinos et al. 2013 ODE's RHS. 
% New version. Avoids 0/0 cases to simplify simulations where species 
% are removed or added in the middle of the run

function dx = Valdovinos2013_rhs(~,x, metadata)

nz_pos  = metadata.nz_pos ;
e       = metadata.e ;
mu_p    = metadata.mu_p ;
mu_a    = metadata.mu_a ;
c       = metadata.c ;
b       = metadata.b ;
u       = metadata.u ;
w       = metadata.w ;
Beta    = metadata.Beta ;
G       = metadata.G ;
g       = metadata.g ;
phi     = metadata.phi ;
tau     = metadata.tau ;
epsilon = metadata.epsilon ;

[p, N, a, Alpha] = unpack(x, metadata) ;
%p(indRemP)=0;% We force extinct species to zero to avoid resucitation due to stiff integrator
%N(indRemP)=0;
%a(indRemA)=0;

% Model's specific computation begins here

sigma = diag(sparse(p.*epsilon)) * Alpha ; %sigma nxm sparse
sigma = sigma * diag(sparse(1./(sum(sigma)+realmin))) ;

Gamma = g .* (1 - u'*p - w.*p + u.*p) ;% non sparse

tmp = Alpha * diag(sparse(a.*tau)) ; %tmp nxm sparse

dp = ( (Gamma .* e .* sum(sigma .* tmp, 2)) - mu_p) .* p;
tmp = (diag(sparse(N)) * tmp) .* b ;
da = sum(c .* tmp, 1)' - mu_a .* a ;
dN = Beta .* p - phi.*N - sum(tmp, 2) ;

% Adaptive dynamics starts here

% Fitness function
DH = diag(sparse(N)) * sparse(c.*b) ; %nxm sparse
DH(Alpha<0)=-DH(Alpha<0) ;

wavg = sum(Alpha.*DH) ; %Weights for average. nxm sparse

%This is the replicator equation
dAlpha = Alpha.*DH - Alpha*diag(sparse(wavg)) ;
dAlpha = dAlpha*diag(sparse(G)) ;


% Now pack the answer
dx = full([dp; dN; da; dAlpha(nz_pos)]) ;
