clear all
close all
clc
load cars
v = [];

%% Annotations

s = whos('-file', 'cars');
for k=1:length(s)
    eval(sprintf('v = [v, annotate(%s)];', s(k).name));
end

%% Centering + scaling

centers = sum(v, 1) / size(v, 1);
vc = v - centers(ones(1, size(v, 1)), :);
ns = sqrt(diag((conj(vc)') * vc)');
vs = vc ./ ns(ones(1, size(v, 1)), :);

%% GPA

S = vs * (conj(vs)');
[V, D] = eigs(S);
m = V(:, 1);
va = (ones(size(v, 1), 1) * ((conj(vs)') * m)') .* vs;

%% PCA

C = cov(va');
[V, D] = eig(C);
[d, i] = sort(diag(D)', 'descend');
V = V(:, i);

s = d * triu(ones(size(D))) / sum(d);
n = find((s >= 0.95), 1);

%% Visualization

for k=1:3
    figure
    drawCarShape(m);
    drawCarShape(m - 3 * sqrt(d(k)) * V(:, k), 2) ;
    drawCarShape(m + 3 * sqrt(d(k)) * V(:, k), 3) ;
end

