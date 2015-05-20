%% DAG clustering experiments with toy DAGs
% Fig. 5 in (Malmi, Tatti, and Gionis, 2015)

% Copyright (c) 2015 Eric Malmi

%%
addpath('extra/')
%%
tic
% Choose either var=1 (vary p_{remove}) or var=2 (vary N_{swaps})
var = 1;
var_strs = {'pr','nis'};

if var == 1
    grid1 = [0:0.05:0.95]; % p_remove
elseif var == 2
    grid1 = [0:1:20];%20]; % n_intra_swaps
end

% Toy DAG densities
p_edges = [0.03 0.08];
p_edge_strs = {'Sparse','Dense'};

k = 5; % Number of clusters
N = 100; % Number of dags
L = 50; % Number of nodes
n_reps = 10; % Number of repetitions

methods = {'median','greedy','hierarchical_complete','seed_DAG_bl'};

res = zeros(n_reps, length(grid1), length(p_edges), length(methods));
res_ari = zeros(n_reps, length(grid1), length(p_edges), length(methods));
% Number of restarts of the k-means with a new random initialization
n_restarts = 10;

tic
fprintf('\n');
for i = 1:length(grid1)
    for e = 1:length(p_edges)
        fprintf('\nval = %.2f, e=%d\n',grid1(i),e);
        n_swaps = 20;
        if var == 1
            nis = 0;
            p_remove = grid1(i);
        elseif var == 2
            nis = grid1(i);
            p_remove = 0;
        end
        p_edge = p_edges(e);
        p_add_new = p_edge*p_remove/(1-p_edge+p_remove*p_edge);
        p = 1/2;
        q = 1/4;

        for r = 1:n_reps
            [As,clusters,plainAs,seedAs] = generate_toy_dags(...
                N, k, L, n_swaps, nis, p_edge, p_remove, p_add_new);

            %drawNetwork('-adjMat',As{1},'-layout',Treelayout);
            %drawNetwork('-adjMat',As{2},'-layout',Treelayout);
            %drawNetwork('-adjMat',As{3},'-layout',Treelayout);
            %drawNetwork('-adjMat',As{4},'-layout',Treelayout);
            %drawNetwork('-adjMat',As{floor(N/k)+1},'-layout',Treelayout);
            %keyboard

            initializations = randi(k,length(As),n_restarts);
            for m = 1:length(methods)
                method = methods{m};

                % Cluster
                if strncmp(method,'hierarchical',12)
                    pdist_mat = zeros(1,N*(N-1)/2);
                    ord = randperm(N);
                    clusters_ord = clusters(ord);
                    for j1 = 1:N
                        for j2 = j1+1:N
                            idx = (j1-1)*(N-j1/2)+j2-j1;
                            pdist_mat(idx) = dag_dist(As{ord(j1)}, ...
                                                      As{ord(j2)}, p, q);
                        end
                    end
                    Z = linkage(pdist_mat,'complete');
                    pred = cluster(Z,'MaxClust',k);
                    ari = rand_index(clusters_ord,pred);

                    res_ari(r,i,e,m) = ari;
                    res(r,i,e,m) = 0;
                    fprintf('%s, ari: %.3f\n',methods{m},ari);
                else
                    min_err = Inf;
                    min_pred = -1;
                    min_n_iter = -1;
                    min_active_clusters = -1;
                    min_centers = -1;
                    for j = 1:n_restarts
                        [pred,err,n_iter,is_removed,centers] = ...
                            graph_k_means(As, k, p, q, 30, ...
                                          initializations(:,j), ...
                                          clusters, method, seedAs, false);
                        if err < min_err
                            min_err = err;
                            min_pred = pred;
                            min_n_iter = n_iter;
                            min_active_clusters = sum(~is_removed);
                            min_centers = centers;
                        end
                    end
                    ari = rand_index(clusters,min_pred);
                    res_ari(r,i,e,m) = ari;
                    res(r,i,e,m) = min_err;
                    avg_center_nodes = 0;
                    for j = 1:k
                        avg_center_nodes = avg_center_nodes + ...
                            sum(sum(min_centers(:,:,j)));
                    end
                    avg_center_nodes = avg_center_nodes / k;
                    fprintf(['%s\terr: %.1f, ari: %.3f, steps: %d, ' ...
                             'clusters: %d, avg_c_n: %.1f\n'], ...
                            methods{m}, min_err, ari, min_n_iter, ...
                            min_active_clusters, avg_center_nodes);
                end
                %keyboard
                %drawNetwork('-adjMat',centers(:,:,1))
            end
        end

        n_edges = 0;
        for j = 1:N
            n_edges = n_edges + sum(As{j}(:));
        end
        fprintf(' %d\n',floor(n_edges/N));
    end
end
fprintf('\n');
toc

%% Plot the results
plot_before(851+10*var,0.8,10,6)
styles = {'b-o', 'r-x', 'k-+','b--*', 'r--s', 'k--v'};
plot_w = length(p_edges);
plot_h = 2;

methods_display = methods;
for i = 1:length(methods_display)
    methods_display{i} = strrep(methods_display{i},'_complete',' ');
    methods_display{i} = strrep(methods_display{i},'_',' ');
end
methods_display{4} = 'optimal';

for e = 1:length(p_edges)
    for m = 1:length(methods)
        if strcmp(methods{m},'hierarchical_average')
            continue
        end
        xs = grid1;
        subplot(plot_h,plot_w,e)
        plot(xs, mean(res_ari(:,:,e,m),1), styles{m}), hold on
        if var == 1
            set(gca,'XTick',[0 0.25 0.5 0.75 1])
            xlabel('p_{remove}')
            title(sprintf('%s DAGs, N_{swaps}=%d',p_edge_strs{e},0))
        elseif var == 2
            xlabel('N_{swaps}')
            title(sprintf('%s DAGs, p_{remove}=%d',p_edge_strs{e},p_remove))
        end
        ylabel('Adjusted Rand Index')
        ylim([-0.05 1.05])
    end

    if e == 2
        subplot(plot_h,plot_w,2)
        legend(methods_display,'Location','SouthWest')
    end
end
