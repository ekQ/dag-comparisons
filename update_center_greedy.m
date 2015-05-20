function [center,cost,edges] = update_center_greedy(graphs,p,q,nV)
%UPDATE_CENTER_GREEDY Computes DAG aggregation for the given graphs.
%
%   [CENTER,COST,EDGES] = UPDATE_CENTER_GREEDY(GRAPHS,P,Q) aggregates
%   graphs given as a (N-by-N-by-K) matrix, where N is the number of edges
%   and K the number of matrices. Returns the aggregated DAG (center), cost
%   (currently not supported, returns -1), and edge list for the center.
%
%   [CENTER,COST,EDGES] = UPDATE_CENTER_GREEDY(GRAPHS,P,Q,NV) aggregates
%   graphs given as a cell array of M-by-2 edge lists. NV is the number of
%   existing vertices.

% Copyright (c) 2015 Eric Malmi

if nargin == 2 || isempty(q)
    q = 0;
end

if nargin == 4 % We have a cell array of edge lists
    edge_counts = zeros(nV, nV);
    all_edges = [];
    M = length(graphs);
    for i = 1:M
        edge_list = graphs{i};
        for j = 1:size(edge_list,1)
            if edge_counts(edge_list(j,1),edge_list(j,2)) == 0
                all_edges = [all_edges; edge_list(j,:)];
            end
            edge_counts(edge_list(j,1),edge_list(j,2)) = ...
                edge_counts(edge_list(j,1),edge_list(j,2)) + 1;
        end
    end
    % Go through unique edges and calculate the regret
    regrets = zeros(size(all_edges,1),1);
    for idx = 1:size(all_edges,1)
        i = all_edges(idx,1);
        j = all_edges(idx,2);
        regrets(idx) = p*(2*edge_counts(i,j) + 2*edge_counts(j,i) - M) ...
            - edge_counts(j,i) + q*(M-edge_counts(i,j)-edge_counts(j,i));
    end
    candidates = all_edges(regrets > 0,:);
    regrets(regrets <= 0) = [];
    [~, order] = sort(regrets,'descend');
    candidates = candidates(order,:);

else % We have adjacency matrices
    graphs_tr = permute(graphs,[2 1 3]);
    ij = sum(graphs,3);
    ji = sum(graphs_tr,3);
    %ij_and_ji_0 = sum(graphs==0 & graphs_tr==0,3);
    %c = p * (ij + ji - ij_and_ji_0) - ji + q * ij_and_ji_0;
    M = size(graphs,3);
    c = p * (2*ij + 2*ji - M) - ji + q * (M - ij - ji);

    [i,j] = find(c>0);
    candidates = [i j];
    values = c(c>0);
    [~, order] = sort(values,'descend');
    candidates = candidates(order,:);
    nV = size(ij,1);
end

center = zeros(nV,nV);
% The transitively closed version of the center
closed_center = zeros(nV,nV);
for l = 1:size(candidates,1)
    i = candidates(l,1);
    j = candidates(l,2);
    if closed_center(j,i) == 0 && center(i,j) == 0
        center(i,j) = 1;
        if closed_center(i,j) == 0
            closed_center(i,j) = 1;
            % Update the transitive closure
            in_nodes = find(closed_center(:,i));
            out_nodes = find(closed_center(j,:));
            closed_center(i,out_nodes) = 1;
            closed_center(in_nodes,j) = 1;
            for ii = in_nodes
                closed_center(ii,out_nodes) = 1;
            end
        end
    end
end
cost = -1;
if ~graphisdag(sparse(center))
    fprintf('Non-cyclic graph!\n')
    keyboard
end

[i,j] = find(center);
edges = [i j];
