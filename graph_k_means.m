function [new_clusters, new_err, n_iterations, is_removed, centers] = ...
    graph_k_means(As, k, p, q, max_iter, initialization, true_clusters, ...
                  method_str, seedAs, debug)
%GRAPH_K_MEANS cluster directed acyclic graphs (DAGs) using a k-means
%   approach.
%
% Input:
%   As              Adjacency matrices cell{input_dag_idx}(node_idx,node_idx)
%   k               Number of clusters
%   p, q            Distance measure parameters
%   max_iter        Maximum number of iterations
%   initialization  Initial cluster assignments vec(dag_idx) = cluster_idx
%   true_clusters   (Optional) True clusters of the matrices. Only needed
%                   for the seed_DAG_bl method.
%   method_str      Which DAG aggregation method (update_center_<methods_str>)
%                   is used. Options: 'greedy' / 'median' / 'random' /
%                   'seed_DAG_bl'.
%   seedAs          (Optional) Adjacency matrices for seed DAGs. Only
%                   needed for the seed_DAG_bl method.
%   debug           Print debug information (true/false).
%
% Output:
%   new_clusters    Cluster indices of the input DAGs
%   new_err         Sum of distances between DAGs and their closest centers
%   n_iterations    How many iterations the algorithm took to converge
%   is_removed      Array indicating which clusters have been removed
%   centers         Cluster centers

% Copyright (c) 2015 Eric Malmi

n_dags = length(As);
% Cluster indices for all input dags
clusters = initialization;
% Number of vertices in a DAG
nV = size(As{1},1);

% Edge lists for the dags
AsE = cell(length(As),1);
for idx = 1:length(As)
    [i,j] = find(full(As{idx}));
    AsE{idx} = [i j];
end

% Edge lists for the seed dags
seedAsE = cell(length(seedAs),1);
for idx = 1:length(seedAs)
    [i,j] = find(full(seedAs{idx}));
    seedAsE{idx} = [i j];
end

% Initial update of the cluster centers
if ~strcmp(method_str,'greedy')
    clusts = cell(k,1); % Dags belonging to different clusters
end
clusts_e = cell(k,1); % Edge lists of dags belonging to different clusters
n_dags_per_cluster = zeros(1,k); % Current number of dags in each cluster
if debug
    fprintf('Starting k-means\n')
end
for i = 1:n_dags
    fullAs = full(As{i});
    if ~graphisdag(sparse(fullAs))
        fprintf('Non-dag input\n')
        keyboard
    end
    idx = n_dags_per_cluster(initialization(i)) + 1;
    if ~strcmp(method_str,'greedy')
        clusts{initialization(i)}(:,:,idx) = fullAs;
    end
    clusts_e{initialization(i)}{idx} = AsE{i};
    n_dags_per_cluster(initialization(i)) = ...
        n_dags_per_cluster(initialization(i)) + 1;
end
% Cluster center adjacency matrices
centers = zeros([size(As{1},1) size(As{1},2) k]);
% Cluster center edge lists
centers_e = cell(k,1);

% Keep track of which clusters are still alive
is_removed = zeros(k,1);
new_err2 = 0;
if debug
    fprintf('Initial assignments\n')
end
for i = 1:k
    if ~isempty(clusts_e{i})
        if strcmp(method_str, 'random')
            [centers(:,:,i),~,centers_e{i}] = update_center_greedy(...
                clusts_e{i}, p, q, nV);
        elseif strcmp(method_str, 'seed_DAG_bl')
            centers(:,:,i) = seedAs{i};
            centers_e{i} = seedAsE{i};
        else
            if ~strcmp(method_str,'greedy')
                [centers(:,:,i),~,centers_e{i}] = feval(...
                    ['update_center_' method_str], clusts{i}, p, q, ...
                    clusts_e{i});
            else
                [centers(:,:,i),~,centers_e{i}] = feval(...
                    ['update_center_' method_str], clusts_e{i}, p, q, nV);
            end
            cost = 0;
            new_err2 = new_err2 + cost;
        end
    else
        %fprintf('Cluster %d removed!\n',l);
        centers(:,:,i) = Inf*ones(size(As{1}));
        centers_e{i} = [];
        is_removed(i) = 1;
        %centers(i,:,:) = As{randi(length(As))};
    end
    if debug
        fprintf('Updated center %d\n',i)
    end
end
%fprintf('Iteration 0, error2 %.2f\n',new_err2)
if debug
    n_dags_per_cluster
    edge_counts = get_edge_counts(centers)
end

% Handle methods 'random' and 'seed_DAG_bl' separately (actually they don't
% belong to this file in the first place!)
if strcmp(method_str, 'random') || strcmp(method_str, 'seed_DAG_bl foo')
    new_err = 0;
    for i = 1:n_dags
        if strcmp(method_str, 'random')
            d = dag_dist(centers(:,:,clusters(i)), full(As{i}), p, q, ...
                         centers_e{clusters(i)}, AsE{i});
        elseif strcmp(method_str, 'seed_DAG_bl')
            d = dag_dist(centers(:,:,true_clusters(i)), full(As{i}), p, ...
                         q, centers_e{true_clusters(i)}, AsE{i});
        end
        new_err = new_err + d;
    end
    n_iterations = 0;
    new_clusters = clusters;
    if strcmp(method_str, 'seed_DAG_bl')
        new_clusters = true_clusters;
    end
    return
end

% The actual k-means
err = Inf;
% How many times the error has increased consequtively
n_incs = 0;
for i = 1:max_iter
    if ~strcmp(method_str,'greedy')
        clusts = cell(k,1);
    end
    clusts_e = cell(k,1);
    n_dags_per_cluster = zeros(1,k); % Current number of dags in each cluster
    new_err = 0;
    new_clusters = zeros(n_dags,1);
    % Assign to clusters
    for j = 1:n_dags
        min_dist = Inf;
        min_idx = -1;
        n_tiers = 1; % Number of clusters within the same distance
        for l = 1:k
            d = dag_dist(centers(:,:,l), full(As{j}), p, q, ...
                         centers_e{l}, AsE{j});
            if d < min_dist
                min_dist = d;
                min_idx = l;
                n_tiers = 1;
            elseif d == min_dist
                % Take this cluster instead of the current min_idx with
                % probability 1/(n_tiers+1)
                n_tiers = n_tiers + 1;
                if rand < 1/n_tiers
                    min_idx = l;
                end
            end
        end
        idx = n_dags_per_cluster(min_idx) + 1;
        if ~strcmp(method_str,'greedy')
            clusts{min_idx}(:,:,idx) = full(As{j});
        end
        clusts_e{min_idx}{idx} = AsE{j};
        n_dags_per_cluster(min_idx) = n_dags_per_cluster(min_idx) + 1;
        new_clusters(j) = min_idx;
        new_err = new_err + min_dist;
    end
    if debug
        n_dags_per_cluster
    end

    if i == 1 && strcmp(method_str, 'seed_DAG_bl')
        break
    end

    %fprintf('Step %d, reconstruction error %f\n',i,new_err);
    if all(new_clusters == clusters)
        if debug
            fprintf('Same clusters\n')
        end
        break
    elseif sum(new_clusters ~= clusters) <= 2
        if debug
            fprintf('Almost same clusters\n')
        end
        break
    elseif new_err == err
        if debug
            fprintf('Clustering cost did not change.\n')
        end
        break
    elseif new_err > err
        n_incs = n_incs + 1;
        if n_incs == 1
            fprintf('  Larger error at %d!\n', i)
            new_clusters = clusters;
            break
        end
    else
        n_incs = 0;
        err = new_err;
        clusters = new_clusters;
        if debug
            fprintf('Iteration %d, error %.2f\n',i,err)
        end
    end
    
    % Update centers
    new_err2 = 0;
    for l = 1:k
        if ~isempty(clusts_e{l})
            if ~strcmp(method_str,'greedy')
                [centers(:,:,l),~,centers_e{l}] = feval(...
                    ['update_center_' method_str], clusts{l}, p, q, ...
                    clusts_e{l});
            else
                [centers(:,:,l),~,centers_e{l}] = feval(...
                    ['update_center_' method_str], clusts_e{l}, p, q, nV);
            end
            cost = 0;
            new_err2 = new_err2 + cost;
        else
            fprintf('    Cluster %d removed!\n',l);
            centers(:,:,l) = Inf*ones(size(As{1}));
            is_removed(l) = 1;
        end
    end
    %fprintf('Iteration %d, error2 %.2f\n',i,new_err2)
    
    if debug
        edge_counts = get_edge_counts(centers)
    end
    if i == max_iter
        fprintf('Maximum number of iterations reached!\n');
    end
end
n_iterations = i;
if debug
    fprintf('--iters %d, err %.2f, clusts %d\n',i,new_err,k-sum(is_removed))
end
%{
% Visualize centers
for i = 1:k
    if ~is_removed(i)
        i
        drawNetwork('-adjMat',squeeze(centers(i,:,:)),'-layout',Treelayout);
    end
end
%}
%fprintf('\n');


function c = get_edge_counts(centers)
c = zeros(1,size(centers,3));
for i = 1:size(centers,3)
    c(i) = sum(sum(centers(:,:,i)));
end

function c = get_member_counts(clusts)
c = zeros(1,length(clusts));
for i = 1:length(clusts)
    c(i) = size(clusts{i},3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist = graph_dist_deltaE(center,graph)
dist = sum(abs(center(:) - graph(:)));
%dist = sum((center(:) - graph(:)).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [center,cost] = update_center_mean(graphs,p)
center = mean(graphs,3);
cost = -1;
