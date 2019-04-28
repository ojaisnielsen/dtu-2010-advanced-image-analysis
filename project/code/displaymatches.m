function [] = displaymatches(matches, varargin)

ims = varargin;
sizecomposite = [0,0];
offset = 0;
for i=1:size(ims,2),
    sizecomposite(2) = sizecomposite(2) + size(ims{i}, 2);
    sizecomposite(1) = max(sizecomposite(1), size(ims{i}, 1));
    offsets(i) = offset;
    offset = offset + size(ims{i}, 2);
end
composite = zeros([sizecomposite, 3]);
for i=1:size(ims,2),
    composite(1:size(ims{i}, 1), offsets(i)+1:offsets(i)+size(ims{i}, 2), :) = ims{i};
end

figure
%colormap('gray');
image(uint8(composite));
hold on

colors = ['r', 'g', 'b'];

for i=1:size(matches,3)-1,
    for j=1:size(matches,2),
        line([offsets(i) + matches(1, j ,i), offsets(i+1) + matches(1, j ,i+1)], [matches(2, j ,i), matches(2, j ,i+1)], 'color', colors(mod(i, size(colors,2))+1));
    end
end

end
