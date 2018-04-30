function [I, J, A] = erdosrenyi(n, p)
% Generate a random Erdos-Renyi graph with n nodes and edge probability p.
%
% [I, J, A] = erdosrenyi(n, p)
% 
% Returns a list of edges (I(k), J(k)) for a random, undirected Erdos-Renyi
% graph with n nodes and edge probability p. A is the adjacency matrix.
%
% I(k) <= J(k) for all k, i.e., all(I<=J) is true.
%
% Nicolas Boumal, UCLouvain, Aug. 7, 2012.

    X = rand(n);
    mask = X <= p;
    X( mask) = 1;
    X(~mask) = 0;
    X = triu(X, 1);

    % A is the adjacency matrix
    A = X + X';
    
    [I, J] = find(X);

end
