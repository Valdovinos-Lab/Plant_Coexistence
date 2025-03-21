function metadata = create_metadata(vectG,In,muAP,i)

% Creates structure with parameter values and initial conditions 
rng(i);
[plant_qty, animal_qty]=size(In);
B=sparse(In);% nxm boolean matrix than describes the network in terms of which
    %    animal can polinize which plant. B(i,j)=1 iff plant i is polinized by
    %    animal j. Try to make this matrix sparse to save memory.
    % (plant_qty=n, animal_qty=m)

nz_pos = find(B); % where the actual links are.

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
c=uniform_rand(mC,vC,plant_qty,animal_qty).*B;
b=uniform_rand(mb,vb,plant_qty,animal_qty).*B;

e=uniform_rand(mE,vE,plant_qty,1);
u=uniform_rand(mU,vU,plant_qty,1);
Beta=uniform_rand(mB,vB,plant_qty,1);
G=uniform_rand(mG,vG,animal_qty,1).*vectG';
g=uniform_rand(mg,vg,plant_qty,1);
mu_a=uniform_rand(mmA,vmA,animal_qty,1);
mu_p=uniform_rand(mmP,vmP,plant_qty,1);
w=uniform_rand(mw,vw,plant_qty,1);
phi=uniform_rand(mphi,vphi,plant_qty,1);
epsilon=uniform_rand(mepsilon,vepsilon,plant_qty,1);

tau=uniform_rand(mtau,vtau,animal_qty,1);


%Give initial state
mz=0.5; vz=1e-1;

yzero=uniform_rand(mz,vz,2*plant_qty+animal_qty,1);

initial_plants=yzero(1:plant_qty);
initial_rewards=yzero(plant_qty+1:2*plant_qty);

initial_animals=yzero(2*plant_qty+1:2*plant_qty+animal_qty);
initial_alphas=B;

%Normalization and packing.
initial_alphas=initial_alphas*diag(sum(initial_alphas).^(-1));
initial_alphas=initial_alphas(nz_pos) ;

% Locate zeros in the Jacobian to speed up the integrator
J_pattern = J_zero_pattern(In) ;

%Checking that all dimensions and values are coherent
assertP ( isequal(size(e), [plant_qty, 1]) ) ;
assertP ( isequal(size(mu_p), [plant_qty, 1]) ) ;
assertP ( isequal(size(mu_a), [animal_qty, 1]) ) ;
assertP ( isequal(size(B), size(c)) ) ;
assertP ( isequal(size(B), size(b)) ) ;
assertP ( isequal(size(u),[plant_qty,1]) ) ;
assertP ( isequal(size(w),[plant_qty,1]) ) ;
assertP ( isequal(size(Beta),[plant_qty,1]) ) ;
assertP ( isequal(size(G),[animal_qty,1]) ) ;
assertP ( isequal(size(g),[plant_qty,1]) ) ;
assertP ( isequal(size(phi),[plant_qty,1]) ) ;

assertP(all(mu_p>0));
assertP(all(mu_a>0));
assertP(all(c>=0));
assertP(all(b>=0));
assertP(all(u>=0));
assertP(all(w>=0));
assertP(all(g>=0));
assertP(all(G>=0));
assertP(all(phi>0));

assertP( all(e) ) ;
assertP( all(mu_p) ) ;
assertP( all(mu_a) ) ;
assertP( all(all(not(c) | B)) ) ;
assertP( all(all(not(b) | B)) ) ;
assertP( all(u) ) ;
assertP( all(w) ) ;
assertP( all(Beta) ) ;
%assertP( all(G) ) ;
assertP( all(g) ) ;
assertP( all(phi) ) ;

metadata = struct( 'plant_qty' , plant_qty,...
                             'animal_qty', animal_qty, ...
                             'nz_pos', nz_pos, ...
                             'e', e, ...
                             'mu_p', mu_p, ...
                             'mu_a', mu_a, ...
                             'c', c, ...
                             'b', b, ...
                             'u', u, ...
                             'w', w, ...
                             'Beta', Beta, ...
                             'G', G, ...
                             'g', g, ...
                             'phi', phi, ...
                             'tau', tau, ...
                             'epsilon', epsilon, ...
                             'B', B, ...
                             'p0', initial_plants, ...
                             'R0', initial_rewards, ... 
                             'a0', initial_animals, ... 
                             'alphas0', initial_alphas,...
                             'J_pattern', J_pattern) ;
end
