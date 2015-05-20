function [As,clusters,plainAs,seedAs] = generate_toy_dags(N,K,L,n_swaps,n_intra_swaps,p_edge,p_remove,p_add_new,dag_order)
%GENERATE_TOY_DAGS Generates directed acyclic graphs
%
% Input:
%   N           Number of dags.
%   K           Number of clusters.
%   L           Total number of items.
%   n_swaps     How many time node indices are randomly swapped for a new
%               cluster (control the similarity of clusters).
%   n_intra_swaps How many swaps within a cluster.
%   p_edge      Probability of creating an edge (i,j) when generating
%               planted clusters.
%   p_remove    Probability of removing an edge when modifying planted
%               clusters.
%   p_add_new   Probability of adding a new edge when modifying planted
%               clusters.
%
% Output:
%   As          Cell array of the transitive closures of the adjancency 
%               matrices of the N dags. Indexing: As{n}(i,j)
%   clusters    N cluster indices (1-K)
%   plainAs     Cell array of the adjancency matrices of the N dags.
%               Indexing: plainAs{n}(i,j)

% Copyright (c) 2015 Eric Malmi

if nargin < 9
    dagord = 1:N;
elseif strcmp(dag_order,'random')
    dagord = randperm(N);
end


As = cell(N,1);
plainAs = cell(N,1);
nPerClust = ceil(N/K);
clusters = zeros(N,1);
seedAs = cell(K,1);
for k = 1:K
    % Modify the edge indices of the planted DAG
    clust_order = 1:L;
    for i = 1:n_swaps
        ord = randperm(L);
        temp = clust_order(ord(1));
        clust_order(ord(1)) = clust_order(ord(2));
        clust_order(ord(2)) = temp;
    end
    
    % Create the planted DAG
    plantA = zeros(L,L);
    if p_edge < 1
        for i = 1:L
            for j = i+1:L
                if rand < p_edge
                    plantA(i,j) = 1;
                end
            end
        end
    else % p_edge equals to the number of edges
        candidates = find(triu(ones(L),1));
        candidates = candidates(randperm(length(candidates)));
        plantA(candidates(1:p_edge)) = 1;
    end
    %{
    % Calculate the transitive closure (tc)    
    tc_plantA = plantA;
    for i = L-1:-1:1
        tc_temp = tc_plantA(i,:);
        incomers = find(tc_plantA(:,i));
        tc_plantA(incomers,:) = tc_plantA(incomers,:) ...
                                | repmat(tc_temp,length(incomers),1);
    end
    %sum(tc_plantA(:)~=tc_plantA2(:))
    %keyboard
    %assert(all(tc_plantA(:)==tc_plantA2(:)));
    %}
    
    % Generate the modifications of the planted clusters
    for j = 1:nPerClust
        dag_idx = (k-1)*nPerClust+j;
        dag_idx = dagord(dag_idx);
        newA = plantA;
        % Remove some existing edges
        rems = rand(L^2,1) < p_remove;
        newA(rems) = 0;
        % Add some new edges
        for i2 = 1:L
            for j2 = i2+1:L
                if rand < p_add_new
                    newA(i2,j2) = 1;
                end
            end
        end
        swapped_order = clust_order;
        for i = 1:n_intra_swaps
            ord = randperm(L);
            temp = swapped_order(ord(1));
            swapped_order(ord(1)) = swapped_order(ord(2));
            swapped_order(ord(2)) = temp;
        end
        newA = newA - diag(diag(newA));
        plainAs{dag_idx} = newA(swapped_order,swapped_order);
        if ~graphisdag(sparse(plainAs{dag_idx}))
            keyboard
        end
        
        % Calculate the transitive closure (tc)
        tc_newA = transitive_closure(newA);
        As{dag_idx} = tc_newA(swapped_order,swapped_order);
        clusters(dag_idx) = k;
        assert(graphisdag(sparse(As{dag_idx})))
    end
    
    % Transitive closure of the seed DAG
    plantA = transitive_closure(plantA);
    seedAs{k} = plantA(clust_order,clust_order);
end
