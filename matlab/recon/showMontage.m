function [f1, f2] = showMontage(img, bkg)
% SHOWMONTAGE displays nD matrix as 2D flattened montage,
%   optionally overlaying on background image
%
%   Usage: showMontage(img, bkg)
%
%       where img is an array of images (if complex, magnitude is used)
%             bkg is a 2D image to use as underlay (optional)
%
%   See also MONTAGE, IMAGESC, VIEWIMAGE
%
%   06/2019, Keith Michel

%%
if ~nargin,      help(mfilename); return,  end
if nargin<2,     bkg = []; end
img = double(img);
bkg = double(bkg);
if ~isreal(img), img = abs(img); end
if ~isreal(bkg), bkg = abs(bkg); end
if ~ismatrix(bkg)
    error('showMontage:bkgSize', ...
        'Input bkg must be 2D but is %d dim array', ndims(bkg))
end

%%
imgSize = size(img);
img = reshape(img, imgSize(1), imgSize(2), []);
nIm = size(img, 3);
nAx = ceil(sqrt(nIm));
nAy = nAx;
if nIm <= (nAx^2)-nAx
    nAx = nAx-1;
end

%%
imgFlat = zeros(nAx*imgSize(1),nAy*imgSize(2));
mask    = zeros(nAx*imgSize(1),nAy*imgSize(2));

idx = 1;
for ii = 1:nAx
    for jj = 1:nAy
        if idx > nIm
            break
        end
        xx = (1:imgSize(1)) + imgSize(1)*(ii-1);
        yy = (1:imgSize(2)) + imgSize(2)*(jj-1);
        imgFlat(xx,yy) = img(:,:,idx);
        mask(xx,yy) = ones(imgSize(1),imgSize(2));
        idx = idx + 1;
    end
end

f1 = figure;
set(gcf, 'position', [0 0 600 600])
set(gca, 'position', [0 0 1 1])
handle = imshow(imgFlat, []);
axis off image
set(handle, 'AlphaData', mask)
if ~isempty(bkg)
    movegui(gcf, 'east')
    f2 = figure;
    set(gcf, 'position', [0 0 600 600])
    set(gca, 'position', [0 0 1 1])
    bkg = repmat(bkg, nAx, nAy);
    overlayImages(bkg, imgFlat, gray(64), hot(64), gca);
    movegui(gcf, 'west')
else
    f2 = [];
    movegui(gcf, 'center')
end

