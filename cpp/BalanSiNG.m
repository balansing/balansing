
function [X] = BalanSiNG(L, E, A, alpha, gamma)
% L: target recursion depth |V| = 2^L
% E: target number of edges
% A: initial adjacency matrix, i.e., A = P + M
% gamma: parameter for noise
% alpha: parameter for weight splitting

sign_method = -1; % stochastic sign decision
rand_seed = 2;

N = 2^L;
S = [1, N];
D = [1, N];

[X] = GenerateEdgesBalanSiNG(L, A, S, D, alpha, gamma, E, rand_seed);

end
