%% Optimize distance measure parameters p and q using Friedman-Rafsky test
% Fig. 3 in (Malmi, Tatti, and Gionis, 2015)

% Copyright (c) 2015 Eric Malmi

%%
addpath('extra/')
%%
vary_p = false; % If false, vary q
% Values for p/q
values = [0.05:0.05:0.95];

% Choose either var=1 (vary p_{remove}) or var=2 (vary N_{swaps})
var = 1;
var_strs = {'pr','nis'};

if var == 1
    grid1 = [0.7]; % p_rems
elseif var == 2
    grid1 = [15]; % n_intra_swaps
end
p_edges = [0.03 0.08];
q_ratios = [0.85 0.45];
p_edge_strs = {'Sparse','Dense'};

N = 100; % Number of dags
L = 50; % Number of nodes
n_reps = 30; % Number of repetitions

% Results
res = zeros(n_reps, length(values), length(grid1), length(p_edges));
tic
fprintf('\n');
for i = 1:length(grid1)
    for e = 1:length(p_edges)
        for qi = 1:length(values)
            if mod(qi,5) == 1
                fprintf('\n%s: val = %.2f, q: %.3f\n', p_edge_strs{e}, ...
                        grid1(i), values(qi));
            end
            if var == 1
                nis = 0;
                p_remove = grid1(i);
            elseif var == 2
                nis = grid1(i);
                p_remove = 0;
            end
            p_edge = p_edges(e);
            % Set edge adding prob so that by expectation the number of
            % edges does not change
            p_add_new = p_edge*p_remove/(1-p_edge+p_remove*p_edge);
            if vary_p
                p = values(qi);
                q = q_ratios(e)*p;
            else
                p = 0.5;
                q = values(qi);
            end
            
            for r = 1:n_reps
                if var == 1
                    n_swaps = 0;
                else
                    n_swaps = 20;
                end
                [As,clusters,plainAs,seedAs] = generate_toy_dags(...
                    N, 2, L, n_swaps, nis, p_edge, p_remove, p_add_new, ...
                    'random');
                % Reference matrix which indicates if two DAGs belong to
                % same distribution 
                refmat = zeros(N);
                D = zeros(N);
                for ii = 1:N
                    for jj = ii+1:N
                        D(ii,jj) = dag_dist(As{ii},As{jj},p,q);
                        if clusters(ii) ~= clusters(jj)
                            refmat(ii,jj) = 1;
                            refmat(jj,ii) = 1;
                        end
                    end
                end
                mst = graphminspantree(sparse(D'));
                edge_idxs = find(mst);
                % Number of edges between DAGs from different clusters
                crossedges = sum(refmat(edge_idxs));
                res(r,qi,i,e) = crossedges;
            end
        end
    end
end
fprintf('\n');
toc

%% Plot the results
figure(802), clf
plot_w = 2;
plot_h = 1;
ymin = [0 10];
ymax = [20 35];
for i = 1:length(grid1)
    for j = 1:length(p_edges)
        ys = mean(res(:,:,1,j),1);
        subplot(plot_h,plot_w,j+(i-1)*plot_w)
        plot(values,ys,'o-')
        ylim([ymin(j) ymax(j)])
        if vary_p
            xlim([0 1])
            xlabel('p')
        else
            xlim([0 0.5])
            xlabel('q')
        end
        ylabel('Test statistic R')
        if var == 1
            if vary_p
                title(sprintf('%s DAGs, q=%.2f*p', p_edge_strs{j}, ...
                              q_ratios(j)))
            else
                title(sprintf('%s DAGs, p=0.5',p_edge_strs{j}))
            end
        else
            title(sprintf('%s DAGs, p_{edge}=%.2f, nis=%d', ...
                          p_edge_strs{j}, p_edges(j), grid1(i)))
        end
    end
end
