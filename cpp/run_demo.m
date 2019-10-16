mexCompile;
A = [0.57, 0.19; 0.19, 0.05];
L = 12; % level
E = 24186; % # of edges to be inserted
alpha = 0.85; % weight parameter % to be changes
gamma = 0.1; % noise parameter

[X] = BalanSiNG(L, E, A, alpha, gamma);

balanced_idx = [0, 3, 5, 6, 8, 10] + 1;
unbalanced_idx = [1, 2, 4, 9, 7, 11] + 1;

pos_idx = find(X(:, 3) > 0);
neg_idx = find(X(:, 3) < 0);

pos_ratio = length(pos_idx)/size(X, 1);
neg_ratio = length(neg_idx)/size(X, 1);

[counts] = SignedDirectedTriangleEnumeration(X);

bal_ratio = sum(counts(balanced_idx)) / sum(counts);
unbal_ratio = sum(counts(unbalanced_idx)) / sum(counts);

fprintf('Ratio of positive edges: %.4f\n', pos_ratio);
fprintf('Ratio of negative edges: %.4f\n', neg_ratio);
fprintf('Ratio of balanced triangles: %.4f\n', bal_ratio');
fprintf('Ratio of unbalanced triangles: %.4f\n', unbal_ratio');
