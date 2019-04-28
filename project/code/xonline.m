function x=xonline(A,B, y)

x = -1;
if (A(2,1) ~= B(2,1)),
    u = (y-A(2,1)) / ((B(2,1)-A(2,1)));
    if (u >= 0 && u <= 1),
        x = A(1,1) + u * (B(1,1)-A(1,1));
    end
end

end