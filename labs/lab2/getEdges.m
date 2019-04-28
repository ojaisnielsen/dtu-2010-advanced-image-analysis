function edges = getEdges(im, t)

r = 3 * round(sqrt(t));
x = [-r:r];

uniDimGauss = exp(-(x.^2)/(2*t))/sqrt(2*pi*t);
uniDimGaussDiffx = -(x/t);
uniDimGaussDiffxx = ((-1/t) + (x/t).^2);
uniDimGaussDiffxxx = ((3*x/(t^2)) - ((x/t).^3));
uniDimGaussDiffy = uniDimGaussDiffx(1 + 2*r:-1:1)';
uniDimGaussDiffyy = uniDimGaussDiffxx(1 + 2*r:-1:1)';
uniDimGaussDiffyyy = uniDimGaussDiffxxx(1 + 2*r:-1:1)';

L = uniDimGauss(ones(1, 2 * r + 1),:) .* uniDimGauss(ones(1, 2 * r + 1),:)';

Lx = uniDimGaussDiffx(ones(1, 2 * r + 1),:) .* L;
LxConv = filter2(Lx, im);
Lxx = uniDimGaussDiffxx(ones(1, 2 * r + 1),:) .* L;
LxxConv = filter2(Lxx, im);
Lxxx = uniDimGaussDiffxxx(ones(1, 2 * r + 1),:) .* L;
LxxxConv = filter2(Lxxx, im);

Ly = uniDimGaussDiffy(:,ones(1, 2 * r + 1)) .* L;
LyConv = filter2(Ly, im);
Lyy = uniDimGaussDiffyy(:,ones(1, 2 * r + 1)) .* L;
LyyConv = filter2(Lyy, im);
Lyyy = uniDimGaussDiffyyy(:,ones(1, 2 * r + 1)) .* L;
LyyyConv = filter2(Lyyy, im);

Lxy = uniDimGaussDiffx(ones(1, 2 * r + 1),:) .* Ly;
LxyConv = filter2(Lxy, im);
Lxxy = uniDimGaussDiffy(:,ones(1, 2 * r + 1)) .* Lxx;
LxxyConv = filter2(Lxxy, im);
Lxyy = uniDimGaussDiffx(ones(1, 2 * r + 1),:) .* Lyy;
LxyyConv = filter2(Lxyy, im);


LvvConv = (LxConv.^2) .* LxxConv + 2 * LxConv .* LyConv .* LxyConv + (LyConv.^2) .* LyyConv;
LvvvConv = (LxConv.^3) .* LxxxConv + 3 * (LxConv.^2) .* LyConv .* LxxyConv + 3 * (LyConv.^2) .* LxConv .* LxyyConv + (LyConv.^3) .* LyyyConv;



LvvP=(LvvConv>0);
Lvv0=xor(LvvP,circshift(LvvP,[0,1])) | xor(LvvP,circshift(LvvP,[1,0]));
edges=Lvv0 & (LvvvConv < 0);

figure
colormap('gray');

imagesc(edges);

% figure
% colormap('gray');
% 
% subplot(4, 4, 1)
% imagesc(L)
% subplot(4, 4, 5)
% imagesc(Lx)
% subplot(4, 4, 6)
% imagesc(Ly)
% subplot(4, 4, 9)
% imagesc(Lxx)
% subplot(4, 4, 10)
% imagesc(Lxy);
% subplot(4, 4, 11)
% imagesc(Lyy);
% subplot(4, 4, 13)
% imagesc(Lxxx);
% subplot(4, 4, 14)
% imagesc(Lxxy);
% subplot(4, 4, 15)
% imagesc(Lxyy);
% subplot(4, 4, 16)
% imagesc(Lyyy);



end