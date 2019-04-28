function result = alphaexp(cost, k, current, x, y, beta)

a = zeros(size(cost, 1),1);
for i=1:size(cost, 1),
    a(i,1) = cost(i, current(i,1));
end
a((current == k)) = inf;
b = cost(:, k);
b((current == k)) = 0;

sourceCost = max(0, b - a);
sinkCost = max(0, a - b);

a = -beta * (current(x) == current(y));
b = -beta * (current(x) == k);
c = -beta * (current(y) == k);
d = -beta * ones(size(x));

sinkCost(y) = sinkCost(y) + c - d;
sourceCost(x) = sourceCost(x) + c - a;
directCost = b + c - a - d;
reverseCost = zeros(size(x));

indices = (1:size(cost, 1))';

TerminalWeights = [indices, sourceCost, sinkCost];
EdgeWeights=[x, y, directCost, reverseCost];
[Cut,Flow]=GraphCutMex(size(cost, 1),TerminalWeights,EdgeWeights);

result = k * ones(size(current));
result(Cut) = current(Cut);


end
