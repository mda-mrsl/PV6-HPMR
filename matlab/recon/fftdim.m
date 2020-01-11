function ksp = fftdim(img, dim)
%FFTDIM computes centered forward discrete Fourier tranform
%
%   Usage: ksp = fftdim(img, dim)
%
%       where img is spatial domain data
%             dim specifies the dimension(s) along which to apply fft 
%               if omitted, every dimension is Fourier transformed
%             ksp is the frequency domain data
%
%   See also SHIFTKSPACE, IFFTDIM
%
%   06/2019, Keith Michel

%% Parse inputs
if ~nargin,      help(mfilename); return; end
if nargin<2,     dim = []; end
if isempty(dim), dim = 1:ndims(img); end
if any(dim > ndims(img))
    error('fftDim:dimensions', ...
        'FT requested along dimension not present in %d-D input array', ...
        ndims(img))
end
sz = size(img);
if any(sz(dim) == 1)
    warning('fftDim:singleton', 'FT applied along singleton dimension')
end

%% Compute forward DFT
ksp = img;
for ii = dim
    ksp = fftshift(fft(ifftshift(ksp, ii), [], ii), ii);
end
