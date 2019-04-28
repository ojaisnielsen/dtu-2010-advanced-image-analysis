%% Load
Im87=freadenvit('thika87');
Im89=freadenvit('thika89');
[m,n,p] = size(Im87);

%% Normalize

Im87n = reshape(Im87, n*m, p);
Im87n = Im87n - ones(m*n, 1) * mean(Im87n);
C = cov(Im87n);
Im87n = Im87n * inv(sqrt(diag(diag(C))));
Im87n = reshape(Im87n, m, n, p);

Im89n = reshape(Im89, n*m, p);
Im89n = Im89n - ones(m*n, 1) * mean(Im89n);
C = cov(Im89n);
Im89n = Im89n * inv(sqrt(diag(diag(C))));
Im89n = reshape(Im89n, m, n, p);

%% Differences

I = Im87n - Im89;
Ih=I(2:end,1:end-1,:);
Iv=I(1:end-1,2:end,:);
Id=(Ih+Iv)/2;
I = I(1:end-1,1:end-1,:);

C = cov(reshape(I, (m-1)*(n-1), p));
Cd = cov(reshape(I-Id, (m-1)*(n-1), p));

[V, D] = eig(Cd, C);
[v, i] = sort(diag(D), 'descend');
V = V(:,i);

Y = reshape(reshape(I, (m-1)*(n-1), p) * V, m-1, n-1, p);
Yh=Y(2:end,1:end-1,:);
Yv=Y(1:end-1,2:end,:);
Yd=(Yh+Yv)/2;
Y = Y(1:end-1,1:end-1,:);

corr(reshape(Y, (m-2)*(n-2), p))
diag(corr(reshape(Y, (m-2)*(n-2), p), reshape(Yd, (m-2)*(n-2), p))) - (-(v/2)+1)


%% Threshold

thres = 0.5;
diff=(abs(Y(:,:,3)) / max(max(abs(Y(:,:,3)))) > thres);


%%
figure
imshowrgb(Im87, [3 2 1], 3);
figure
imshowrgb(Im89, [3 2 1], 3);
figure
imagesc(diff);

%% MRF

[H, W] = size(Y(:,:,3));
v = reshape(abs(Y(:,:,3)), H*W, 1);

cost(:,2) = max(v) - v;
cost(:,1) = v;

%%

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


%%
figure
imshowrgb(Im87, [3 2 1], 3);
figure
imshowrgb(Im89, [3 2 1], 3);
figure
imagesc(result);


