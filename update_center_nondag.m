function [center,cost] = update_center_nondag(graphs,p,q)
if nargin == 2
    q = 0;
end
graphs_tr = permute(graphs,[2 1 3]);
ij = sum(graphs,3);
ji = sum(graphs_tr,3);
ij_and_ji_0 = sum(graphs==0 & graphs_tr==0,3);
c = p * (ij + ji - ij_and_ji_0) - ji + q * ij_and_ji_0;
majority = ij>ji | transpose(ij<ji);
majority(c<=0) = 0;
%center = transitive_closure(majority);
center = majority;
cost = -1;