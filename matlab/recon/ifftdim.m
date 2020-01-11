function img = ifftdim(ksp, dim)
%IFFTDIM computes centered inverse discrete Fourier tranform
%
%   Usage: img = ifftdim(ksp, dim)
%
%       where ksp is the frequency domain data
%             dim specifies the dimension(s) along which to apply ifft 
%               if omitted, every dimension is Fourier transformed
%             img is spatial domain data
%
%   See also SHIFTKSPACE, FFTDIM
%
%   06/2019, Keith Michel

%% Parse inputs
if ~nargin,      help(mfilename); return; end
if nargin<2,     dim = []; end
if isempty(dim), dim = 1:ndims(ksp); end
if any(dim > ndims(ksp))
    error('ifftDim:dimensions', ...
        'FT requested along dimension not present in %d-D input array', ...
        ndims(ksp))
end
sz = size(ksp);
if any(sz(dim) == 1)
    warning('ifftDim:singleton', 'FT applied along singleton dimension')
end

%% Compute forward DFT
img = ksp;
for ii = dim
    img = fftshift(ifft(ifftshift(img, ii), [], ii), ii);
end
