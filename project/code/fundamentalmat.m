function F = fundamentalmat(points1, points2, norm)




if norm,
    [points1, Tn1] = normalize(points1);
    [points2, Tn2] = normalize(points2);
end

M = [points2(1,:)'.*points1(1,:)', points2(1,:)'.*points1(2,:)', points2(1,:)', ...
    points2(2,:)'.*points1(1,:)', points2(2,:)'.*points1(2,:)', points2(2,:)', ...
    points1(1,:)', points1(2,:)', ones(size(points1,2),1) ];

[U,S,V] = svd(M,0);
F = reshape(V(:,9),3,3)';

[U,S,V] = svd(F);
F = U*diag([S(1,1), S(2,2), 0])*V';

if norm,
    F = Tn2' * F * Tn1;
end


end

