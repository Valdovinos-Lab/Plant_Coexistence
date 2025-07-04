%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modification 08/03/2019, Ann Arbor
% Cleaning up my codes
% Only run the dynamics, without species invasions or removals
% Runs the dynamics for only one matrix.
% Outputs the whole time-series for any variable as well as final values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [plantsf, Rf, animalsf, alphasf]=IntegrateValdovinos2013_dataset(metadata)

tmax=2000000;
%EUp=2e-2;
%EUa=1e-3;


% Combining all initial variables
initial_state=full([metadata.p0;metadata.R0;metadata.a0;metadata.alphas0]);
tspan = [0 tmax];

%% Integrating the dynamic model
options = odeset('JPattern', metadata.J_pattern,'NonNegative',1:2*metadata.plant_qty+metadata.animal_qty) ;

rhs_with_metadata = @(t, y) Valdovinos2013_rhs(t, y, metadata);

[~, y]=ode15s(rhs_with_metadata, tspan, initial_state, options) ;

yf = y(end,:)';

% Retriving final densities and foraging efforts
[plantsf, Rf, animalsf, alphasf] = unpack(yf, metadata);

plantsf=full(plantsf);
Rf=full(Rf);
animalsf=full(animalsf);
alphasf=full(alphasf);

end
