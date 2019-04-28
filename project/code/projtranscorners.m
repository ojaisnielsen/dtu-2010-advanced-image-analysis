function corners = projtranscorners(x0, y0, x1, y1, T)

corners = [x0, x0, x1, x1; y1, y0, y0, y1; 1, 1, 1, 1];
corners = T * corners;
corners = corners(1:2,:) ./ (ones(2, 1) * corners(3,:));

end