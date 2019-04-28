initvlfeat
 %% Load data

photos = dir('input');
im = [];
K = [];
n = 1;
s=0.5;
inch = 25.4;
for i=1:size(photos),
    parts = regexp(photos(i).name, '\.', 'split');
    if (strcmpi(parts(end), 'jpg')),
        filename = strcat('input/', photos(i).name)
        im{n} = imresize(imread(filename), s);
        info = imfinfo(filename);
        f = info.DigitalCamera.FocalLength;
        cx = size(im{n}, 2) / 2;
        cy = size(im{n}, 1) / 2;
        K(:,:,n) = [s*f*info.XResolution /inch,0,cx;0,s*f*info.YResolution/inch,cy;0,0,1]
        n = n + 1;
    end
end

ims = [2,3];
pair = 1;
i1 = ims(pair);
i2 = ims(pair+1);

%% crop

for i=ims,
    
    figure
    imagesc(im{i});

    [x1,y1] = ginput(1);
    [x2,y2] = ginput(1);
    im{i} = im{i}(min(y1, y2):max(y1,y2), min(x1, x2):max(x1,x2),:);
    K(1,3,i) = K(1,3,i) - min(x1, x2);
    K(2,3,i) = K(2,3,i) - min(y1, y2);
    figure
    imagesc(im{i});
end



%% SIFT features extraction and matches
clear feat desc points matches
for i=1:size(ims,2),
    [feat{i}, desc{i}] = VL_SIFT(single(rgb2gray(im{ims(i)})));
    points{i} = [feat{i}(1:2, :); ones(1, size(feat{i}, 2))]; 
end

for i=1:size(ims,2)-1,
    matches{i} = VL_UBCMATCH(desc{i}, desc{i+1});  
end

%% Harris corners extraction and matches

sigma = 3;
thresh = 100;
radius = 3*sigma;
clear corels points matches
for i=1:size(ims,2),
    [cim, r, c] = harris(rgb2gray(im{ims(i)}), sigma, thresh, radius, 1);
    points{i} = [c'; r'; ones(1, size(r, 1))]; 
end
radius = 4;
for i=1:size(points,2)-1,
    corels = zeros(size(points{i}, 2), size(points{i+1}, 2));
    for p=1:size(points{i}, 2),
        for q=1:size(points{i+1}, 2),
            yp = points{i}(2,p);
            xp = points{i}(1,p);
            yq = points{i+1}(2,q);
            xq = points{i+1}(1,q);
            if not(yp-radius <= 0 || yp+radius > size(im{ims(i)}, 1) || xp-radius <= 0 ||xp+radius > size(im{ims(i)}, 2) || yq-radius <= 0 || yq+radius > size(im{ims(i+1)}, 1) || xq-radius <= 0 ||xq+radius > size(im{ims(i+1)}, 2)),
                corels(p,q) = normcrosscorel(im{ims(i)}(yp-radius:yp+radius, xp-radius:xp+radius, :), im{ims(i+1)}(yq-radius:yq+radius, xq-radius:xq+radius, :));
            end
        end
    end
    [C, best2] = max(corels,[], 1);
    [C, best1] = max(corels,[], 2);
    n=0;
    for p=1:size(points{i}, 2),
        if (best2(best1(p)) == p && corels(p, best1(p)) > 0.5),
            n = n+1;
            matches{i}(:,n) = [p; best1(p)];
        end
    end
end


%% Robust matches and estimation of F

homogconstraint = false;

npoints = points;
if not(homogconstraint),
    for i=1:size(ims,2),   
        [npoints{i}, Tn{i}] = normalize(points{i});
    end    
end

[fullmatches, F, H] = robustmatches(0.99, 0.01, homogconstraint, npoints, matches);

clear nmatches
for i=1:size(ims,2),
    nmatches(:,:,i) = points{i}(:,fullmatches(i,:));
end

if not(homogconstraint),
    for i=1:size(ims,2)-1,
        F(:,:,i) = Tn{i+1}' * F(:,:,i) * Tn{i};
    end
end

for i=1:size(ims,2)-1,
    F(:,:,i) = fundamentalmat(nmatches(:,:,i), nmatches(:,:,i+1), true);
end

displaymatches(nmatches, im{ims});


%% Manual matches

sizecomposite = [0,0];
offset = 0;
for i=ims,
    sizecomposite(2) = sizecomposite(2) + size(im{i}, 2);
    sizecomposite(1) = max(sizecomposite(1), size(im{i}, 1));
    offsets(i) = offset;
    offset = offset + size(im{i}, 2);
end

composite = zeros([sizecomposite, 3]);
for i=ims,
    composite(1:size(im{i}, 1), offsets(i)+1:offsets(i)+size(im{i}, 2), :) = im{i};
end

figure
image(uint8(composite));
hold on

colors = ['r', 'g', 'b'];


select = true;
nmatches = [];
j = 0;
while select,
    j = j+1;
    for i=1:size(ims, 2),
        [x,y] = ginput(1);
        nmatches(:, j, i) = [x - offsets(ims(i)); y; 1]; 
        plot(x,y, 'r+')
    end
    [x,y,button] = ginput(1);
    select = (button == 1);
end   

for i=1:size(ims,2)-1,
    F(:,:,i) = fundamentalmat(nmatches(:,:,i), nmatches(:,:,i+1), true);
end


displaymatches(nmatches, im{ims});




%% Check epipolar lines

p=1

close all
P = nmatches(:,p,1);
Q = nmatches(:,p,2);
imagesc(im{ims(1)});
hold on
plot(P(1), P(2), 'g+');
drawline(F(:,:,1)'* Q, 1, size(im{ims(1)}, 2), 1, size(im{ims(1)}, 1));
figure
imagesc(im{ims(2)});
hold on
plot(Q(1), Q(2), 'g+');
drawline(F(:,:,1) * P, 1, size(im{ims(2)}, 2), 1, size(im{ims(2)}, 1));

p=p+1;

%% Camera matrices estimation from calibrated cameras

K1 = K(:,:,i1);
K2 = K(:,:,i2);

R1 = eye(3);
t1 = zeros(3,1);

E = K2'*F(:,:,pair)*K1;
[U,S,V] = svd(E);
W = [0,-1,0;1,0,0;0,0,1];
t = U(:,3);

R2s(:,:,1) = U * W * V';
R2s(:,:,2) = U * W' * V';
R2s(:,:,3) = U * W * V';
R2s(:,:,4) = U * W' * V';

if (det(R2s(:,:,1)) < 0 && det(R2s(:,:,2)) <0 && det(R2s(:,:,3)) < 0 && det(R2s(:,:,4)) <0),
    [U,S,V] = svd(-E);
    W = [0,-1,0;1,0,0;0,0,1];
    t = U(:,3);

    R2s(:,:,1) = U * W * V';
    R2s(:,:,2) = U * W' * V';
    R2s(:,:,3) = U * W * V';
    R2s(:,:,4) = U * W' * V';
end


t2s(:,:,1) = t;
t2s(:,:,2) = t;
t2s(:,:,3) = -t;
t2s(:,:,4) = -t;

for i=1:4,
    points1 = triangulate(nmatches(:,:,pair), nmatches(:,:,pair+1), K1*[R1,t1], K2 * [R2s(:,:,i), t2s(:,:,i)]);
    depth1(1,i) = mean(sign(points1(3,:)));
    
%     figure
%     plot3(points1(1,:)', points1(2,:)', points1(3,:)', 'o', 'MarkerSize',3,'MarkerEdgeColor','r','MarkerFaceColor','r');
    
    points2 = (R2s(:,:,i)' * points1(1:3,:)) + t2s(:,:,i) * ones(1, size(points1,2));
    depth2(1,i) = mean(sign(points2(3,:)));
    
%     figure
%     plot3(points2(1,:)', points2(2,:)', points2(3,:)', 'o', 'MarkerSize',3,'MarkerEdgeColor','r','MarkerFaceColor','r');    

end
depth1
depth2

i = find((depth1 > 0) .* (depth2 > 0))
R2 = R2s(:,:,i)
t2 = t2s(:,:,i)


%% Camera matrices estimation from 3 views

perm = randperm(size(nmatches, 2));
i3 = setdiff(ims, [i1, i2]);
i3 = i3(1);
[P,X] = vgg_PX_from_6pts_3img(nmatches(:,perm(1:6),[i1, i2, i3]));

s = 3;
[K1, R1, t1] = vgg_KR_from_P(P(:,:,1,s));
[K2, R2, t2] = vgg_KR_from_P(P(:,:,2,s));
t1 = -R1*t1;
t2 = -R2*t2;


%% Rectification


c1 = -R1' * t1;
c2 = -R2' * t2;
xaxis = c2-c1;
yaxis = cross(R1(3,:)', xaxis);
zaxis = cross(xaxis, yaxis);

Krect = (K1 + K1) / 2;
Krect(1,2) = 0;
R = [xaxis'/norm(xaxis); yaxis'/norm(yaxis); zaxis'/norm(zaxis)];
t1rect = - R * c1;
t2rect = - R * c2;

T1 = Krect * R * R1' * inv(K1);
T2 = Krect * R * R2' * inv(K2);

corners1 = projtranscorners(1, 1, size(im{i1}, 2), size(im{i1}, 1), T1);
corners2 = projtranscorners(1, 1, size(im{i2}, 2), size(im{i2}, 1), T2); 
imrectsize1 = [max(corners1(2,:)) - min(corners1(2,:)) + 1, max(corners1(1,:)) - min(corners1(1,:)) + 1];
imrectsize2 = [max(corners2(2,:)) - min(corners2(2,:)) + 1, max(corners2(1,:)) - min(corners2(1,:)) + 1];

C1 = [K1(1,3); K1(2,3); 1];
C1rect = T1 * C1;
C1rect = C1rect / C1rect(3,1);
C2 = [K2(1,3); K2(2,3); 1];
C2rect = T2 * C2;
C2rect = C2rect / C2rect(3,1);

dx1 = C1rect(1,1) - C1(1,1);
dx2 = C2rect(1,1) - C2(1,1);
dy = C1rect(2,1) - C1(2,1);

K1rect = Krect - [0,0,dx1;0,0,dy;0,0,0];
K2rect = Krect - [0,0,dx2;0,0,dy;0,0,0];


s = max([size(im{i1}), size(im{i2})]) / max([imrectsize1, imrectsize2]);
S = [s,0,0;0,s,0;0,0,1];
K1rect = S * K1rect;
K2rect = S * K2rect;

T1 = K1rect * R * R1' * inv(K1);
T2 = K2rect * R * R2' * inv(K2);

corners1 = projtranscorners(1, 1, size(im{i1}, 2), size(im{i1}, 1), T1);
corners2 = projtranscorners(1, 1, size(im{i2}, 2), size(im{i2}, 1), T2); 
boundingbox1 = [min([corners1(1,:), corners2(1,:)]), max([corners1(1,:), corners2(1,:)]); min([corners1(2,:), corners2(2,:)]), max([corners1(2,:), corners2(2,:)])];
boundingbox2 = boundingbox1

im1rect = imtransform(im{i1}, maketform('projective', T1'), 'YData', boundingbox1(2,:), 'XData', boundingbox1(1,:));
im2rect = imtransform(im{i2}, maketform('projective', T2'), 'YData', boundingbox2(2,:), 'XData', boundingbox2(1,:));

corners1 = corners1 - (boundingbox1(:,1) - ones(2,1)) * ones(1, 4);
corners2 = corners2 - (boundingbox2(:,1) - ones(2,1)) * ones(1, 4);

figure
colormap('gray');
imagesc(im1rect);
hold on
line([corners1(1,:), corners1(1,1)], [corners1(2,:), corners1(2,1)], 'color', 'r');

figure
colormap('gray');
imagesc(im2rect);
hold on
line([corners2(1,:), corners2(1,1)], [corners2(2,:), corners2(2,1)], 'color', 'r');

oldcorners1 = corners1;

%% Crop


figure
imagesc(im1rect);

corners1 = []
for i=1:4,
    [x,y] = ginput(1);
    corners1(1, i) = x;
    corners1(2, i) = y;
    line(corners1(1,:), corners1(2,:), 'color', 'r');
end
line([corners1(1,:), corners1(1,1)], [corners1(2,:), corners1(2,1)], 'color', 'r');





%% Rectification error estimation


% im1rect = imread('tsukuba/scene1.row3.col2.ppm');
% im2rect = imread('tsukuba/scene1.row3.col1.ppm');
% corners1 = [1, 1, size(im1rect, 2), size(im1rect, 2); size(im1rect, 1), 1, 1, size(im1rect, 1)];
% corners2 = [1, 1, size(im2rect, 2), size(im2rect, 2); size(im2rect, 1), 1, 1, size(im2rect, 1)];

homogconstraint = false;

clear feat desc points Tn matches tuples
[feat{1}, desc{1}] = VL_SIFT(single(rgb2gray(im1rect)));
points{1} = [feat{1}(1:2, :); ones(1, size(feat{1}, 2))]; 
[feat{2}, desc{2}] = VL_SIFT(single(rgb2gray(im2rect)));
points{2} = [feat{2}(1:2, :); ones(1, size(feat{2}, 2))]; 
    
npoints = points;
for i=1:2,
    if not(homogconstraint),
        [npoints{i}, Tn{i}] = normalize(points{i});
    end
end
matches{1} = VL_UBCMATCH(desc{1}, desc{2});  

[fullmatches, F, H] = robustmatches(0.99, 0.001, homogconstraint, npoints, matches);

clear matches
for i=1:2,
    matches(:,:,i) = points{i}(:,fullmatches(i,:));
end

displaymatches(matches, im1rect, im2rect);

meandisp = round(mean(matches(1,:,2) - matches(1,:,1)));

[XI, YI]= meshgrid(1:size(im1rect,2),1:size(im1rect,1));
X = [matches(1,:,1), oldcorners1(1,:)];
Y = [matches(2,:,1), oldcorners1(2,:)];
Z = [matches(2,:,1) - matches(2,:,2), 0,0,0,0];
ZI = round(griddata(X, Y, Z, XI, YI));
ZI(isnan(ZI)) = 0;

figure
colormap('gray');
imagesc(ZI)

%% Disparity estimation

disparities = meandisp-20:meandisp+20
ndisps = size(disparities, 2);

x0s = zeros(size(im1rect,1), ndisps);
x1s = zeros(size(im1rect,1), ndisps);
y0s = zeros(ndisps,1);
y1s = zeros(ndisps,1);

for y=1:size(im1rect,1),
    x1 = -1*ones(4,1);
    x2 = -1*ones(4,1);
    for i=1:4,
        x1(i,1) = xonline(corners1(:,i), corners1(:,mod(i,4)+1), y);
        x2(i,1) = xonline(corners2(:,i), corners2(:,mod(i,4)+1), y);
    end
    x1 = x1(logical((x1 > 0)));
    x2 = x2(logical((x2 > 0)));
    if (size(x1, 1) > 0 && size(x2, 1) > 0),
        di=0;
        for d=disparities,
            di = di+1;
            x2d = x2 - d;
            if (size(x2d, 1) > 0),
                x0s(y, di) = max(round(max(min(x1), min(x2d))),1);
                x1s(y, di) = min(round(min(max(x1), max(x2d))), size(im1rect,2));
            end
        end
    end
end
for d=1:ndisps,  
    i = find(x0s(:,d) > 0);
    y0s(d,1) = i(1,1);
    y1s(d,1) = i(end,1);
end


radius = 4;

clear dispmapcost
dispmapcost = cell(ndisps);
means = zeros(ndisps,1);

time = cputime;
di = 0;
%parfor d=disparities,
for d=disparities,
% d=disparities(1,1);
    di = di+1
    m = 0;
    n = 0;
    y0 = y0s(di,1)+radius;
    y1 = y1s(di,1)-radius;
    dispmapcost{di} = nan(size(im1rect, 1), size(im1rect, 2));
    for y=y0:y1,
        x0 = x0s(y, di)+radius;
        x1 = x1s(y, di)-radius;
        for x=x0:x1,
            dispmapcost{di}(y,x) = -log(normcrosscorel(im1rect(y-radius:y+radius, x-radius:x+radius, :), im2rect(y-radius-ZI(y,x):y+radius-ZI(y,x), x-radius+d:x+radius+d, :)));
                m = m + dispmapcost{di}(y,x);
                n = n+1;
            
            %dispmapcost{di}(y,x) = exp(-normcrosscorel(im1rect(y-radius:y+radius, x-radius:x+radius), im2rect(y-radius:y+radius, x-radius+d:x+radius+d)));
             %dispmapcost{di}(y,x) = exp(-normcrosscorel(im1rect(y-radius:y+radius, x-radius:x+radius), im2rect(y-ZI(y,x)-radius:y-ZI(y,x)+radius, x-radius+d:x+radius+d)));
             %dispmapcost{di}(y,x) = (im1rect(y,x) - im2rect(y-ZI(y,x), x + d))^2;

        end
    end
    means(di,1) = m/n;
end
cputime-time
% save 'save2.mat' -v7.3

%%
dispmapcostfinal = zeros([size(im1rect, 1), size(im1rect, 2), ndisps]);



unionmask = false(size(im1rect, 1), size(im1rect, 2));
for d=1:ndisps,
    mask = not(isnan(dispmapcost{d}));
    unionmask = not(not(unionmask).*not(mask));
    dispmapcostfinal(:,:,d) = mask .* dispmapcost{d} + means(d,1) * not(mask);
%     mask = (dispmapcostfinal(:,:,d) > log(2));
%     dispmapcostfinal(:,:,d) = means(d,1) * mask + not(mask) .* dispmapcostfinal(:,:,d);
%     dispmapcostfinal(:,:,d) = bitmax * mask + not(mask) .* dispmapcostfinal(:,:,d);
end

[C,I] = min(dispmapcostfinal,[], 3);
figure
colormap('gray');
imagesc(I);


%%


H = size(im1rect, 1);
W = size(im1rect, 2);


vectmask = reshape(unionmask, W*H, 1);

cost = reshape(dispmapcostfinal, W*H, ndisps);
cost = cost(vectmask, :);
vect = reshape(I, W*H, 1);
vect = vect(vectmask);

beta = 50;

indices = (1:H*W)';
maskindices = zeros(size(vect));
maskindices(vectmask) = (1:sum(vectmask))';

% Vertical cliques
xV = indices(mod(indices, H) ~= 0);
xV = xV(logical(vectmask(xV) .* vectmask(xV + 1)));
yV = xV + 1;

xV = maskindices(xV);
yV = maskindices(yV);


% Horizontal cliques
xH = indices((indices + H) <= W*H);
xH = xH(logical(vectmask(xH) .* vectmask(xH + H)));
yH = xH + H;

xH = maskindices(xH);
yH = maskindices(yH);

% Solution
x = [xV; xH];
y = [yV; yH];

for i=1:10
    i
    for k=1:size(cost, 2)                
        vect = alphaexp(cost, k, vect, x, y, beta);
    end
end

result = zeros(H*W, 1);
result(vectmask) = vect;
result = reshape(result, H, W);
figure
colormap('gray')
imagesc(result)

%%


[x, y] = meshgrid(1:size(im1rect,2), 1:size(im1rect,1));
x = reshape(x, 1, size(im1rect,1)*size(im1rect,2));
y = reshape(y, 1, size(im1rect,1)*size(im1rect,2));
% x2 = x + disparities(:, reshape(result, 1, size(im1rect,1)*size(im1rect,2)));
c = reshape(im1rect, 1, size(im1rect,1)*size(im1rect,2), 3);

x = x(vectmask);
y = y(vectmask);

disps = reshape(result, 1, size(im1rect,1)*size(im1rect,2));
x2 = x + (disparities(disps(vectmask)));
c = double(c(1, vectmask, :)) / 255;

points = triangulate([x; y; ones(1, size(x, 2))], ...
    [x2; y; ones(1, size(x, 2))], ...
    K1rect * [R, t1rect], K2rect * [R, t2rect]);


X = points(1,:);
Y = points(2,:);
Z = points(3,:);


mX = min(X);
MX = max(X);
mY = min(Y);
MY = max(Y);
[XI, YI]= meshgrid(mX:(MX - mX +1)/500:MX,mY:(MY - mY +1) /500:MY);

r = c(1,:, 1);
g = c(1,:, 2);
b = c(1,:, 3);
ZI = griddata(X, Y, Z, XI, YI);
RI = griddata(X, Y, r, XI, YI);
GI = griddata(X, Y, g, XI, YI);
BI = griddata(X, Y, b, XI, YI);
mask = not(isnan(ZI));
% ZI(isnan(ZI)) = max(max(ZI(not(isnan(ZI)))));
RI(isnan(RI)) = 0;
GI(isnan(GI)) = 0;
BI(isnan(BI)) = 0;
clear C
C(:,:,1) = RI;
C(:,:,2) = GI;
C(:,:,3) = BI;
% figure
% colormap('gray')
% surf(XI, YI, ZI)
% vrml(gcf,'output.wrl')

figure
% surf(XI, YI, ZI)
surf(XI, YI, ZI, C, 'FaceColor','interp','edgecolor','none','FaceLighting','phong')

% surf(XI, YI, ZI)
% vrml(gcf,'output.wrl')
camlight right

%%







