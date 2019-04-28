function imshowrgb(image, rgb, varargin)

%
% imshowrgb(image, rgb, varargin)
%
% image - 3-D matrix
% rgb   - order of R, G and B, e.g. [4 2 1]
%         defaults to [1 2 3]
% nstd  - strech between mean -/+ nstd stddevs of image bands,
%         defaults to 2

if ndims(image)~=3, error('input image must be 3-D'); end

nstd = 2;
if nargin==3
    nstd = varargin{1};
end
[nr nc nf] = size(image);
image = image(:,:,rgb);
nf = 3;
image = reshape(image,nc*nr,nf);
image = image-repmat(mean(image),nr*nc,1);
stdvar = std(image);
image = reshape(image,nr,nc,nf);
r = image(:,:,1)/(2*nstd*stdvar(1))+0.5; r(r<0)=0; r(r>1)=1;
g = image(:,:,2)/(2*nstd*stdvar(2))+0.5; g(g<0)=0; g(g>1)=1;
b = image(:,:,3)/(2*nstd*stdvar(3))+0.5; b(b<0)=0; b(b>1)=1;
imshow(cat(3,r,g,b))
