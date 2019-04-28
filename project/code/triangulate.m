function X = triangulate(x1, x2, P1, P2)


X = zeros(size(x1, 1), size(x1,2));
for i=1:size(x1,2),
    M = [crossprodmat(x1(:, i)) * P1; crossprodmat(x2(:, i)) * P2];
    [U,S,V] = svd(M,0);
    x = V(:,end);    
    X(:, i) = x(1:end-1, 1) / x(end,1);     
end




end
