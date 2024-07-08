function [D_i, w_degree, degree]=effective_degree(alphasf)

% Calculates effective degree based on information theory metrics of
% biodiversity by Hill 1973:
% D_i^q=[sum_j(alpha_ij/sum_k(alpha_ik))^q]^1/(1-q)

% Bofore: [D_i, D_Vi, w_degree, w_degreeVi, degree]=effective_degree(alphasf, Visits);

% Effective degree w/alpha
q=2;
inside_sum=(alphasf./sum(alphasf,2)).^q;
D_i=sum(inside_sum,2).^(1/(1-q));

% % Effective degree w/Visits (not useful because it has the confounding
% abundances, plus they drive the Di to NaN when plants go extinct)
% sumVisit=diag(sum(Visits,2))*ones(size(Visits));
% inside_sum=(Visits./sumVisit).^q;
% D_Vi=sum(inside_sum,2).^(1/(1-q));

% Weigthed degree w/alpha
w_degree=sum(alphasf,2);

% % Weighted degree w/Visits( not useful because it has the confounding
% abundances, plus they drive the Di to 0 when plants go extinct)
% sumVisit=diag(sum(Visits,2))*ones(size(Visits));
% w_degreeVi=sum(Visits,2);

% degree
degree = sum(alphasf > 0, 2);

end