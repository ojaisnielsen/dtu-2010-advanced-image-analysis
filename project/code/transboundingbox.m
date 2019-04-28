function [x0, y0, x1, y1] = transboundingbox(x0, y0, x1, y1, T)

M = T * [x0, x1, x1, x0; y0, y0, y1, y1; 1,1,1,1];
M = M ./ (ones(3,1) * M(3,:));

x0 = min(M(1,:));
y0 = min(M(2,:));
x1 = max(M(1,:));
y1 = max(M(2,:));

end
