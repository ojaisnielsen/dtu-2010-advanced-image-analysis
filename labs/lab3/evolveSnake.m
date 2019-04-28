function snake = evolveSnake(im, snake, weighting, deviation, alpha, beta)

n = size(snake, 1);
[H, W] = size(im);
gauss = fspecial('gaussian', 3 * round(deviation),deviation);
gaussConv = filter2(gauss, im);
[Px, Py] = gradient(gaussConv);
P = -weighting * (Px.^2 + Py.^2);
[Fx, Fy] = gradient(P);
Fx = -Fx;
Fy = -Fy;

d0 = -2 * alpha - 6 * beta;
d1 = alpha + 4 * beta;
d2 = -beta;


A = diag(d1 * ones(n - 1, 1), 1) + diag(d2 * ones(n - 2, 1), 2);
A(1,n) = d1;
A(1,n-1) = d2;
A(2,n) = d2;
A = A + A' + diag(d0 * ones(n, 1));
M = eye(n) - 1 * A;
% M = sparse(m);
Minv = inv(M);

isConverged = false;
i = 0;
while not(isConverged)
    XY1d=round(snake(:,2))+size(Fx,1)*(round(snake(:,1))-1);
    Fext=[Fx(XY1d) Fy(XY1d)];
    
    snake = Minv * (snake + 1 * Fext);
    snake = [max(min(snake(:,1), W), 1), max(min(snake(:,2), H), 1)];
    if i == 10000
        isConverged = true;
    end
    i = i + 1;
end

end
