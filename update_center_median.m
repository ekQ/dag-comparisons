function [center,cost,edges] = update_center_median(graphs,p,q,edge_lists)
%UPDATE_CENTER_MEDIAN Computes DAG aggregation for the given graphs by
%   selecting the graph with the smallest distance to all others.
%
%   Input: 
%       graphs(start_node_idx, end_node_idx, graph_idx)
%       p, q are distance measure parameters
%       edge_lists (optional) is a cell array of M-by-2 edge lists.

% Copyright (c) 2015 Eric Malmi

if nargin == 2
    q = 0;
end
n = size(graphs,3);
D = zeros(n,n);
for i = 1:n
    for j = i:n
        if nargin > 3
            D(i,j) = dag_dist(graphs(:,:,i),graphs(:,:,j),p,q, ...
                              edge_lists{i},edge_lists{j});
        else
            D(i,j) = dag_dist(graphs(:,:,i),graphs(:,:,j),p,q);
        end
    end
end
D = D + D' - diag(diag(D));
[cost,c_idx] = min(sum(D));
try
    center = graphs(:,:,c_idx);
catch e
    keyboard
end
if ~graphisdag(sparse(center))
    keyboard
end

[i,j] = find(center);
edges = [i j];
