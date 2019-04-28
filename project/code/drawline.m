function [] = drawline(l, x0, x1, y0, y1)
n = 0;
X = [];
Y = [];
if l(1) ~= 0,
    if (-(l(3) + l(1) * x0) / l(2) >= y0 && -(l(3) + l(1) * x0) / l(2) <= y1)
        X = [x0];
        Y = [-(l(3) + l(1) * x0) / l(2)];
        n = n + 1;
    end
    if (-(l(3) + l(1) * x1) / l(2) >= 1 && -(l(3) + l(1) * x1) / l(2) <= y1)
        X = [X, x1];
        Y = [Y, -(l(3) + l(1) * x1) / l(2)];
        n = n + 1;
    end
end
if l(2) ~= 0,
    if (n < 2 && -(l(3)+l(2)*y0)/l(1) >= x0 && -(l(3)+l(2)*y0)/l(1) <= x1)
        X = [X, -(l(3)+l(2)*y0)/l(1)];
        Y = [Y, y0];
        n = n + 1;        
    end
    if (n < 2 && -(l(3)+l(2)*y1)/l(1) >= x0 && -(l(3)+l(2)*y1)/l(1) <= x1)
        X = [X, -(l(3)+l(2)*y1)/l(1)];
        Y = [Y, y1]; 
        n = n+1;
    end    
end
    
line(X, Y, 'color', 'r');



end
