function c = normcrosscorel(x, y)
    nx = single(reshape(x, size(x,1) * size(x,2) * size(x,3), 1));
    ny = single(reshape(y, size(y,1) * size(y,2) * size(y,3), 1));
    nx = nx - mean(nx);
    ny = ny - mean(ny);
    nx = nx / norm(nx);
    ny = ny / norm(ny);
    c = dot(nx, ny);
    c = max(eps, c);
end
