%% DAG Aggregation experiments with toy DAGs
% Fig. 4 in (Malmi, Tatti, and Gionis, 2015)

% Copyright (c) 2015 Eric Malmi

%%
addpath('extra/')
%%
vars = [1 2];
var_strs = {'pr','nis'};

grids = {};
grids{1} = 0:0.05:0.95;
grids{2} = 0:1:20;
grid_lengths = zeros(length(grids),1);
for i = 1:length(grids)
    grid_lengths(i) = length(grids{i});
end
longest_grid = max(grid_lengths);

p_edges = [0.03 0.08];
p_edge_strs = {'Sparse','Dense'};

N = 100; % Number of dags
L = 50; % Number of nodes
n_reps = 30; % Number of repetitions

methods = {'median','greedy','empty_baseline','seed_baseline'};

% Store the results (the matrix might contain extra space due to unequal
% grid lengths, should use cell array instead)
res = zeros(n_reps, longest_grid, length(p_edges), length(vars), ...
             length(methods));

tic
fprintf('\n');
for vi = 1:length(vars)
    var = vars(vi);
    for i = 1:length(grids{vi})
        for e = 1:length(p_edges)
            fprintf('\n%s: val = %.2f\n',p_edge_strs{e},grids{vi}(i));
            if var == 1
                nis = 0;
                p_remove = grids{vi}(i);%0;
            elseif var == 2
                nis = grids{vi}(i);
                p_remove = 0;
            end
            p_edge = p_edges(e);
            p_add_new = p_edge*p_remove/(1-p_edge+p_remove*p_edge);
            p = 1/2;
            q = 1/4;

            for r = 1:n_reps
                [As,~,plainAs,seedAs] = generate_toy_dags(...
                    N,1,L,0,nis,p_edge,p_remove,p_add_new);
                As2 = zeros(L,L,length(As));
                for j = 1:length(As)
                    As2(:,:,j) = As{j};
                end
                for m = 1:length(methods)
                    method = methods{m};
                    if strcmp(method,'seed_baseline')
                        center = seedAs{1};
                    else
                        % Aggregate the DAGs
                        center = feval(['update_center_' method],As2,p,q);
                    end
                    % The q value used for comparing the methods
                    q_comparison = q;
                    sumK = 0;
                    for j = 1:length(As)
                        sumK = sumK + dag_dist(As{j}, center, p, ...
                                               q_comparison);
                    end
                    res(r,i,e,vi,m) = dag_dist(center, seedAs{1}, p, ...
                                               q_comparison);
                    center_edge_n = sum(center(:));
                    if r == 1
                        fprintf(['%s\terr: %.1f, dist2seed: %.1f, ' ...
                                 'center_edges: %d\n'], ...
                                 methods{m},sumK,res(r,i,m),center_edge_n);
                    end
                end
            end
        end
    end

    n_edges = 0;
    for j = 1:N
        n_edges = n_edges + sum(As{j}(:));
    end
    fprintf(' %d\n',floor(n_edges/N));
end
fprintf('\n');
toc

%% Plot the results
figure(702), clf
styles = {'b-o', 'r-x', 'k-+', 'b--*', 'r--s', 'k--v','b-.d', 'r-.^', ...
          'k-.p'};
plot_w = 2;
plot_h = 2;

methods_display = methods;
for i = 1:length(methods_display)
    methods_display{i} = strrep(methods_display{i},'_',' ');
end
methods_display{4} = 'optimal';

for e = 1:length(p_edges)
    for vi = 1:length(vars)
        for m = 1:length(methods)
            xs = grids{vi};
            var = vars(vi);

            plot_idx = (vi-1)*plot_w+e;
            subplot(plot_h,plot_w,plot_idx)
            plot(xs, mean(res(:,1:length(xs),e,vi,m),1), styles{m})
            hold on
            if var == 1
                set(gca,'XTick',[0 0.25 0.5 0.75 1])
                xlabel('p_{remove}')
                title(sprintf('%s DAGs, N_{swaps}=%d',p_edge_strs{e},20))
            elseif var == 2
                xlabel('N_{swaps}')
                title(sprintf('%s DAGs, p_{remove}=%d',p_edge_strs{e},0))
            end
            if e == 1
                ylabel('K(center,seed)')
                ylim([280 330])
            end
        end
    end
end
subplot(plot_h,plot_w,4)
legend(methods_display,'Orientation','horizontal');
