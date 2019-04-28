function v = annotate(im)
%ANNOTATE Summary of this function goes here
%   Detailed explanation goes here
figure
colormap('gray');
imagesc(im)
u = MarkShape();
v = u(:,1)+i*u(:,2);

end
