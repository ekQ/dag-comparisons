%% Cluster users based pairwise artist preferences and predict preferences
% Sec. 7.3 in (Malmi, Tatti, and Gionis, 2015)

% Copyright (c) 2015 Eric Malmi

addpath('extra/')
% If you want to visualize graphs, download graphViz4Matlab
% http://www.mathworks.com/matlabcentral/fileexchange/21652-graphviz4matlab
% and add it to path
%addpath('extra/graphViz4Matlab/')

%% Cyclicity of the preference graphs
prefs = load('music_preferences/music_preferences.dat');
n_all_users = max(prefs(:,1));
artists = importdata('music_preferences/list_of_artists.txt');
n_nodes = 20;
is_dag = [];

n_train_users = 480;
D = {}; % Train dags
D2 = {}; % Test dags
D2plain = {}; % Test dags without the transitive closure
% User indices for the preference dags
user_idxs = [];
user_idxs2 = [];
for i = 1:n_all_users
    links = prefs(prefs(:,1)==i,2:3);
    graph = sparse(links(:,1), links(:,2), ones(size(links,1),1), ...
                   n_nodes, n_nodes);
    is_dag(i) = graphisdag(graph);
    if ~is_dag(i)
        continue
    end
    
    if length(D) < n_train_users
        D{end+1} = transitive_closure(full(graph));
        user_idxs = [user_idxs; i];
    else
        D2{end+1} = transitive_closure(full(graph));
        D2plain{end+1} = graph;
        user_idxs2 = [user_idxs2; i];
    end
end
fprintf('Fraction of cyclic graphs: %.4f\n', sum(~is_dag)/length(is_dag));

n_users = length(D);
n_test_users = length(D2);

%% Plot artist popularities
figure(1), clf
likes = histc(prefs(:,3),1:20);
bar(likes)
set(gca, 'XTick', 1:20, 'XTickLabel',artists)
rotateXLabels(gca, 45)

%% Cluster the train users
%ks = [1:10 15:5:30];
ks = [3];

p = 0.5;
q = 0.48;
methods = {'median','greedy'};

n_restarts = 30;

tic
res = zeros(length(methods),length(ks));
preds = zeros(length(methods),n_users,length(ks));

res_all = zeros(length(methods),n_restarts);
times = zeros(length(methods), n_restarts);
iters = zeros(length(methods), n_restarts);

for ki = 1:length(ks)
    k = ks(ki);
    for m = 1:length(methods)
        method = methods{m};
        % Cluster
        if strncmp(method,'hierarchical',12)
            n_dags = n_users;
            pdist_mat = zeros(1,n_dags*(n_dags-1)/2);
            ord = randperm(n_dags);
            for j1 = 1:n_dags
                for j2 = j1+1:n_dags
                    idx = (j1-1)*(n_dags-j1/2)+j2-j1;
                    pdist_mat(idx) = dag_dist(D{ord(j1)},D{ord(j2)},p,q);
                end
            end
            Z = linkage(pdist_mat,'complete');
            min_pred = cluster(Z,'MaxClust',k);
        else
            min_err = Inf;
            min_pred = -1;
            min_n_iter = -1;
            min_active_clusters = -1;
            min_centers = -1;
            for j = 1:n_restarts
                initialization = randi(k,n_users,1);
                t_beg = tic;
                [pred,err,n_iter,is_removed,centers] = graph_k_means(...
                    D, k, p, q, 30, initialization, [], method, [], false);
                times(m,j) = toc(t_beg);
                res_all(m,j) = err;
                iters(m,j) = n_iter;
                if err < min_err
                    min_err = err;
                    min_pred = pred;
                    min_n_iter = n_iter;
                    min_active_clusters = sum(~is_removed);
                    min_centers = centers;
                end
            end
            preds(m,:,ki) = min_pred';
            avg_center_nodes = 0;
            for j = 1:k
                avg_center_nodes = avg_center_nodes + ...
                    sum(sum(min_centers(:,:,j)));
                % NB: Uncomment the following lines if you wish to
                % visualize the preference DAGs (and pick your k)
                %drawNetwork('-adjMat', min_centers(:,:,j), ...
                %    '-nodeLabels', artists);%, '-layout',Treelayout)
                %keyboard
            end
            avg_center_nodes = avg_center_nodes / k;
            fprintf(['%s\terr: %.1f, steps: %d, clusters: %d, ' ...
                     'avg_c_n: %.1f\n'], methods{m}, min_err, ...
                     min_n_iter, min_active_clusters, avg_center_nodes);
            res(m,ki) = min_err;
        end
    end
end
toc

%% Plot the results
figure(26), clf
plot(ks,res,'-o')
xlabel('# of clusters')
ylabel('Total distance')

%% Predict preferences for test users
% First prepare the cluster matrices for prediction
n_nodes = size(D{1},1);
cmat = zeros(n_nodes, n_nodes, k);
for i = 1:k
    members = find(min_pred==i);
    for j = 1:length(members)
        cmat(:,:,i) = cmat(:,:,i) + D{members(j)};
    end
end

train_edges = 1:13;
n_reps = 30;

all_res = zeros(n_reps, length(train_edges));
all_res_bl = zeros(n_reps, length(train_edges));
for ei = 1:length(train_edges)
    % Number of edges used for finding the best matching cluster
    ne = train_edges(ei);
    for ri = 1:n_reps
        res = [];
        res_bl = [];
        for i = 1:n_test_users
            edges = find(D2plain{i});
            if length(edges) <= ne
                continue
            end
            edges = edges(randperm(length(edges)));
            graph_train = zeros(n_nodes, n_nodes);
            % Add the ne first edges and take the transitive closure
            graph_train(edges(1:ne)) = 1;
            graph_train = transitive_closure(graph_train);
            % And the remaining edges to the test graph
            graph_test = zeros(n_nodes, n_nodes);
            graph_test(edges(ne+1:end)) = 1;
            % However, remove those test edges that are already present in
            % the train graph
            graph_test = graph_test - (graph_test & graph_train);

            % Calc distances to different clusters
            dists = zeros(k,1);
            for j = 1:k
                dists(j) = dag_dist(graph_train,min_centers(:,:,j),p,q);
            end
            pred_mat = zeros(n_nodes, n_nodes);
            bl_mat = zeros(n_nodes, n_nodes);
            for j = 1:k
                % Pick the nearest cluster
                pred_mat = pred_mat + (dists(j)==min(dists))*cmat(:,:,j);
                %pred_mat = pred_mat + dists(j)*cmat(:,:,j);
                bl_mat = bl_mat + cmat(:,:,j);
            end

            % Find the pairs we want to predict
            [xs,ys] = find(graph_test);
            for j = 1:length(xs)
                if pred_mat(xs(j),ys(j)) > pred_mat(ys(j),xs(j))
                    res = [res; 1];
                elseif pred_mat(xs(j),ys(j)) == pred_mat(ys(j),xs(j))
                    res = [res; 0.5];
                else
                    res = [res; 0];
                end
            end
            for j = 1:length(xs)
                if bl_mat(xs(j),ys(j)) > bl_mat(ys(j),xs(j))
                    res_bl = [res_bl; 1];
                elseif bl_mat(xs(j),ys(j)) == bl_mat(ys(j),xs(j))
                    res_bl = [res_bl; 0.5];
                else
                    res_bl = [res_bl; 0];
                end
            end
        end
        all_res(ri,ei) = mean(res);
        all_res_bl(ri,ei) = mean(res_bl);
    end
    fprintf('Edges: %d. accuracy: %.4f, baseline: %.4f, n: %d\n', ne, ...
            mean(all_res(:,ei)), mean(all_res_bl(:,ei)), length(res));
end

%% Plot results
figure(54), clf
plot(train_edges, mean(all_res_bl), 'k-x'), hold on
plot(train_edges, mean(all_res), '-o'), hold on
xlabel('# of train preferences')
xlim([1 13])
set(gca, 'XTick', [1:2:13])
ylabel('Accuracy')
legend('Majority vote (baseline)', 'Majority vote within cluster', ...
       'Location', 'SouthEast')
