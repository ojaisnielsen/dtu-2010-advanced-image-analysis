function V = nullvector(M)

% [V, D] = eig(M);
% V = V(:, 1);

[U,S,V] = svd(M,0);
V = V(:,end);

end