function g = transitive_closure(graph)
%TRANSITIVE_CLOSURE computes the transitive closure of the given adjacency
%   matrix.
g = graph;
N = size(graph,1);
g = g - diag(diag(g)) + diag(ones(N,1));
g = g>0;
g = (g^N)>0;
g = g - diag(diag(g));