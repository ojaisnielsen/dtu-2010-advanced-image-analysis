%% Load images

nIm = 100;
Images=[];
for cIm=0:(nIm-1),
    if(cIm<10)
        name=sprintf('ukbench0000%d.jpg',cIm);
    elseif(cIm<100)
        name=sprintf('ukbench000%d.jpg',cIm);
    end
    Images{cIm+1}=rgb2gray(imread(name));
end
%% Find keypoints
SIFTdescr=[];
for cIm=1:nIm,
    imwrite(Images{cIm},'HaaScript.pgm');
    [image, descrips, locs] = sift('HaaScript.pgm');
    SIFTdescr{cIm}=descrips;
end
%%
features = []
for cIm=1:nIm,
    features = [features; SIFTdescr{cIm}];
end

%% K-means
k=50;

[centers,mincenter,mindist,q2,quality] = kmeans2(features,k);

%% Image descriptors
imageDescriptors = [];
for cIm=1:nIm,
    imageDescriptors = [imageDescriptors; zeros(1, k)];
    for cSift=1:size(SIFTdescr{cIm}, 1),
        [m, argm] = min(sqrt(sum(((ones(k,1) * SIFTdescr{cIm}(cSift, :)) - centers).^2,  2)));
        imageDescriptors(cIm, argm) = imageDescriptors(cIm, argm) + 1;
    end
end

%% Test

testIm = 1;

X = ones(nIm - 1, 1) * imageDescriptors(testIm, :);
Y = imageDescriptors((1:nIm ~= testIm), :);
S = X+Y;
D = (X - Y).^2;
D((S == 0)) = 0;
S((S == 0)) = 1;
dist = sum(D ./ S, 2);
plot(dist)






