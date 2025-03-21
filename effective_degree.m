function [D_j, wd_i]=effective_degree(Alpha)

% Effective degree of pollinators w/alpha
% Calculates effective degree based on information theory metrics of
% biodiversity by Hill 1973:
% D_i^q=[sum_j(alpha_ij/sum_k(alpha_ik))^q]^1/(1-q)
q=2;
inside_sum=(Alpha./sum(Alpha)).^q;
D_j=sum(inside_sum).^(1/(1-q));

% Weigthed plant degree w/alpha
wd_i=sum(Alpha,2);


end