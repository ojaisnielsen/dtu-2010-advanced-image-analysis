function H = homography(P, Q, norm)

if norm,
    [P, Tp] = normalize(P);
    [Q, Tq] = normalize(Q);
end

M = [];
for i=1:size(P,2),
    M = [M; kron(crossprodmat(Q(:, i)), P(:,i)')];
end

[U,S,V] = svd(M,0);
H = reshape(V(:,end), 3, 3)';

if norm,
    H = inv(Tq) * H * Tp;
end

end