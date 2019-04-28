nFrames=10; Draw=1;
[x,y,P,Lines]=Film(nFrames,Draw);

%% 8-point algo

cFrame1 = 1;
cFrame2 = 5;

B = [x(cFrame2,:)' .* x(cFrame1,:)', x(cFrame2,:)' .* y(cFrame1,:)', x(cFrame2,:)', y(cFrame2,:)' .* x(cFrame1,:)', y(cFrame2,:)' .* y(cFrame1,:)', y(cFrame2,:)', x(cFrame1,:)', y(cFrame1,:)', ones(size(x(cFrame2,:)'))];  

[U,S,V] = svd(B,0);
F = reshape(V(:,9),3,3)';

l2 = F*[x(cFrame1,1); y(cFrame1,1); 1];
l1 = F'*[x(cFrame2,1); y(cFrame2,1); 1];


%% Draw epipolar line on first frame

drawline(l1, min(x(cFrame1,:)), max(x(cFrame1,:)), min(y(cFrame1,:)), max(y(cFrame1,:)));

%% Draw epipolar line on second frame

drawline(l2, min(x(cFrame2,:)), max(x(cFrame2,:)), min(y(cFrame2,:)), max(y(cFrame2,:)));

%% Ransac
prob=0.99;
data = RanLine(10,15);
thres = 1;

data = [data; ones(1, size(data, 2))];
plot(data(1,:), data(2,:), 's')
hold on

converged = false;
maxScore=0;
best=[];
n = 0;
N = inf; 
while N > n,
    perm = randperm(size(data, 2));
    m = ones(3);
    m(1:3,1:2) = data(:,perm(1:2));
    l = (m'\[0;0;1])';
    inliers = (abs(l*data)/norm(l) <= thres);
    if (sum(inliers) > maxScore),
        maxScore = sum(inliers);
        best = data(:, inliers);
    end
    
    eps = 1-(sum(inliers)/size(data, 2));       
    N = log(1-prob)/log(1 - (1 - eps)^2);    
    n = n+1;    
end

l = polyfit(best(1,:), best(2,:), 1);
x=[floor(min(data(1,:))), ceil(max(data(1,:)))];
plot(x,polyval(l,x))

%% 3d inference

image = imread('petergade.png');
imagesc(image)

Q = [0, 6.7+1.98+3.96+0.76, 0, 6.7+1.98+3.96+0.76;
    0, 0, 6.1, 6.1];

Q = [100 * Q; ones(1, size(Q,2))];


[Px,Py] = ginput(size(Q,2));
P = [Px'; Py'; ones(1, size(Q,2))];

%%

[Pn,Tp] = normalize(P);
[Qn,Tq] = normalize(Q);


M = [];
for i=1:size(P,2),
    M = [M; [zeros(1,3), -Qn(3,i)*(Pn(:,i)'), Qn(2,i)*(Pn(:,i)')]; [Qn(3,i)*(Pn(:,i)'), zeros(1,3), -Qn(1,i)*(Pn(:,i)')]; [-Qn(2,i)*(Pn(:,i)'), Qn(1,i)*(Pn(:,i)'), zeros(1,3)]];
end

[U,S,V] = svd(M,0);
H = reshape(V(:,9),3,3);
H = Tp' * H * inv(Tq');

figure
Tr=maketform('projective', H);
WarpIm=imtransform(image,Tr,'YData',[0 Q(2,3)],'XData',[0 Q(1,2)]);
imagesc(WarpIm)

[x,y] = ginput(2);
