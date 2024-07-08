%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modification 08/03/2019, Ann Arbor
% Cleaning up my codes
% Only run the dynamics, without species invasions or removals
% Runs the dynamics for only one matrix.
% Outputs the whole time-series for any variable as well as final values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [plantsf, nectarf, animalsf, alphasf,network_metadata]=IntegrateValdovinos2013_dataset(vectG,In,muAP)

global J_pattern network_metadata
tmax=10000;
%EUp=2e-2;
%EUa=1e-3;

[m, n]=size(In);
B=sparse(In);

% Parameters of the uniform distribution from where the parameters of the
% dynamic model are drawn:
varp=1e-1;% variance of plant parameters
vara=0;% variance of animal parameters
mC=0.2; vC=vara;
mE=0.8; vE=varp;
mb=0.4; vb=vara;
mU=0.06; vU=varp;
mw=1.2; vw=varp;
mB=0.2; vB=varp;
mG=2; vG=vara;
mg=0.4; vg=varp;
mphi=0.04; vphi=varp;
mtau=1; vtau=vara;
mepsilon=1; vepsilon=0;% Keeping this variance 0, as we are considering epsilon 1
                       % for all coexistence calculations (i.e., ignoring it)
vmA=vara; vmP=varp;

if muAP==1
    mmA=0.05; mmP=0.001; % high pollinator mortality 
elseif muAP==2
    mmA=0.001; mmP=0.02; % high plant mortality
elseif muAP==3
    mmA=0.001; mmP=0.001; % low plant and animal mortality
elseif muAP==4
    mmA=0.03; mmP=0.005; % high plant and animal mortality
end

% Parameters are drawn from uniform distribution (see Valdovinos et al.
% 2018, Nature Communications for complete description of the model and
% parameter definition)

% (10%meanP)-meanP+(10%meanP); (0.01%meanA)-meanA+(0.01%meanA)
c=uniform_rand(mC,vC,m,n).*B;
b=uniform_rand(mb,vb,m,n).*B;

e=uniform_rand(mE,vE,m,1);
u=uniform_rand(mU,vU,m,1);
Beta=uniform_rand(mB,vB,m,1);
G=uniform_rand(mG,vG,n,1).*vectG';
g=uniform_rand(mg,vg,m,1);
mu_a=uniform_rand(mmA,vmA,n,1);
mu_p=uniform_rand(mmP,vmP,m,1);
w=uniform_rand(mw,vw,m,1);
phi=uniform_rand(mphi,vphi,m,1);
epsilon=uniform_rand(mepsilon,vepsilon,m,1);

tau=uniform_rand(mtau,vtau,n,1);

%Create structure 
network_metadata = create_metadata(B, e, mu_p, mu_a, c, b, u, w, Beta, G, g, phi, tau, epsilon) ;

%Give initial state
mz=0.5; vz=1e-1;

yzero=uniform_rand(mz,vz,2*m+n,1);

initial_plants=yzero(1:m);
initial_nectar=yzero(m+1:2*m);

initial_animals=yzero(2*m+1:2*m+n);
initial_alphas=B;

%Normalization and packing.
initial_alphas=initial_alphas*diag(sum(initial_alphas).^(-1));

initial_alphas=initial_alphas(network_metadata.nz_pos) ;

% Combining all initial variables
initial_state=full([initial_plants;initial_nectar;initial_animals;initial_alphas]);
tspan = [0 tmax];

%% Integrating the dynamic model
options = odeset('JPattern', J_pattern,'NonNegative',1:2*m+n) ;
[t, y]=ode15s(@Valdovinos2013_rhs,tspan,initial_state, options) ;

yf = y(end,:)';

% Retriving final densities and foraging efforts
[plantsf, nectarf, animalsf, alphasf] = unpack(yf, network_metadata);

plantsf=full(plantsf);
nectarf=full(nectarf);
animalsf=full(animalsf);
alphasf=full(alphasf);

end
