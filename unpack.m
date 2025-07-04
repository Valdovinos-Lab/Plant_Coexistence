function [plants, R, animals, alphas] = unpack(y, metadata)

n = metadata.plant_qty ;
m = metadata.animal_qty ;

plants=y(1:n,1) ;
R=y(n+1:2*n,1);
animals=y(2*n+1: 2*n+m,1);

alphas=sparse( n, m ) ;
alphas(metadata.nz_pos) = y(2*n+m+1:end,1) ;

end
