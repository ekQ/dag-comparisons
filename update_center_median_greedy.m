function [center,cost] = update_center_median_greedy(graphs,p,q)
%UPDATE_CENTER_MEDIAN_GREEDY Computes DAG aggregation for the given graphs
%   by selecting the graph with the smallest distance to all others and
%   then greedily adding new edges if it improves the center.

% Copyright (c) 2015 Eric Malmi

if nargin == 3
    q = 0;
end
graphs_tr = permute(graphs,[2 1 3]);
ij = sum(graphs,3);
ji = sum(graphs_tr,3);
ij_and_ji_0 = sum(graphs==0 & graphs_tr==0,3);
c = p * (ij + ji - ij_and_ji_0) - ji + q * ij_and_ji_0;

candidates = find(c(:)>0);
values = c(candidates);
[values, order] = sort(values,'descend');
candidates = candidates(order);

% Find the closest and set it as initial center
n = size(graphs,3);
D = zeros(n,n);
for i = 1:n
    for j = i+1:n
        D(i,j) = dag_dist(graphs(:,:,i),graphs(:,:,j),p,q);
    end
end
D = D + D';
[cost,c_idx] = min(sum(D));
center = graphs(:,:,c_idx);

% Improve the center
for l = 1:length(candidates)
    [i,j] = ind2sub(size(ij),candidates(l));
    if center(j,i) == 0 && center(i,j) == 0
        center(i,j) = 1;
        % Update the transitive closure
        center(i,find(center(j,:))) = 1;
        center(find(center(:,i)),j) = 1;
    end
end
cost = -1;
if ~graphisdag(sparse(center))
    fprintf('Non-cyclic center (closest greedy)\n');
    %keyboard
end
