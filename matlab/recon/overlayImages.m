function [axs, ims] = overlayImages(undr, over, cmapu, cmapo, ax)
% OVERLAYIMAGES overlays images with transparency and color
%
%   Usage: [axs, ims] = overlayImages(undr, over, [cmapu, cmapo, ax])
%
%       where undr is a 2D image to underlay (if complex, magnitude is used)
%       	  over is a 2D image to overlay (if complex, magnitude is used)
%       	  cmapu is the colormap to use for underlay, default gray(256)
%       	  cmapo is the colormap to use for overlay, default hot(256)
%       	  ax is a handle to an axes on which to draw images
% 
%       	  axs contains the handles of underlay and overlay axes
%       	  ims contains the handles of underlay and overlay images
%
%   See also MONTAGE, IMAGESC, VIEWIMAGE, SHOWMONTAGE, VIEWOVER
%
%   06/2019, Keith Michel

%%
if nargin<2,      help(mfilename); return,  end
if nargin<3,      cmapu = []; end
if nargin<4,      cmapo = []; end
if nargin<5,      ax = []; end
if ~isreal(undr), undr = abs(undr); end
if ~isreal(over), over = abs(over); end
if ~ismatrix(undr) || ~ismatrix(over)
    error('overlayImages:dimensions', 'Images must be 2-dimensional')
end
if isempty(cmapu),  cmapu = gray(256); end
if isempty(cmapo),  cmapo = hot(256); end
if isempty(ax)    
    figure; 
    ax = axes('position', [0 0 1 1]);
else
    cla(ax);
end

%%
undr = nscale(undr);
axes(ax)
if verLessThan('matlab', '8.4') % Guessing on version number here...
    undr = 1 + round(undr * (size(cmapu,1)-1));
    im1 = subimage(undr, cmapu);
else
    im1 = imagesc(undr);
    colormap(ax, cmapu)
end
axis off image
ax2 = axes;
im2 = imagesc(over);
colormap(gca, cmapo)
axis off image
set(im2, 'AlphaData', nscale(over))
axs = [ax, ax2];
ims = [im1, im2];
setappdata(ax, 'position_listener', linkprop(axs,'Position'))
