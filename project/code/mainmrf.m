close all
clear all
clc
load brain
colormap('gray')
imagesc(brain)
muCSF = 32;
muGM = 47;
sdCSF = 5;
sdGM = 4;
muWM = 55;
sdWM = 3;

[H, W] = size(brain);
v = reshape(brain, H*W, 1);

cost = [(log(sdCSF^2) + ((v - muCSF).^2) / (sdCSF^2)) / 2, (log(sdGM^2) + ((v - muGM).^2) / (sdGM^2)) / 2, (log(sdWM^2) + ((v - muWM).^2) / (sdWM^2)) / 2];

%% 1-clique
sourceCost = max(0, cost(:,2) - cost(:,1));
sinkCost = max(0, cost(:,1) - cost(:,2));

TerminalWeights = [(1:H*W)', sourceCost, sinkCost];

EdgeWeights=[1, 2, 0, 0];

[Cut,Flow]=GraphCutMex(H*W,TerminalWeights,EdgeWeights);

result = ones(H*W, 1);
result(Cut', :) = 0;
result = reshape(result, H, W);
figure
colormap('gray')
imagesc(result)

%% 2-clique

beta = 3;

b = cost(:,2);
a = cost(:,1);

sourceCost = max(0, b - a);
sinkCost = max(0, a - b);

indices = (1:H*W)';

% Vertical cliques
xV = indices(mod(indices, H) ~= 0);
yV = xV + 1;

% Horizontal cliques
xH = indices((indices + H) <= W*H);
yH = xH + H;

% Solution
x = [xV; xH];
y = [yV; yH];

a = -beta * ones(size(x));
b = zeros(size(x));
c = zeros(size(x));
d = -beta * ones(size(x));

sinkCost(y) = sinkCost(y) + c - d;
sourceCost(x) = sourceCost(x) + c - a;
directCost = b + c - a - d;
reverseCost = zeros(size(x));

TerminalWeights = [indices, sourceCost, sinkCost];
EdgeWeights=[x, y, directCost, reverseCost];
[Cut,Flow]=GraphCutMex(H*W,TerminalWeights,EdgeWeights);

result = ones(H*W, 1);
result(Cut', :) = 0;
result = reshape(result, H, W);
figure
colormap('gray')
imagesc(result)

%% Alpha expansion

beta = 3;

indices = (1:H*W)';

% Vertical cliques
xV = indices(mod(indices, H) ~= 0);
yV = xV + 1;

% Horizontal cliques
xH = indices((indices + H) <= W*H);
yH = xH + H;

% Solution
x = [xV; xH];
y = [yV; yH];

result = init(cost);
for i=1:10
    for k=1:size(cost, 2)
        result = alphaexp(cost, k, result, x, y, beta);
    end
end

result = reshape(result, H, W);
figure
colormap('gray')
imagesc(result)


