function imgOut = calcMzImg(imgIn, noi, exc)
%CALCMZIMG rescales MR magnitude image data to values proportional to Mz
%
%   Usage: imgOut = calcMzImg(imgIn, noi, exc)
%
%       where imgIn is Fourier transformed magnitude image data
%               4th dimension corresponds to excitation angles in exc
%             noi contains magnitude image noise (must be FT'ed)
%             exc is a vector of excitation angles in degrees
%               if omitted, default of 90 degrees is used
%             imgOut is the rescaled magnitude image data
%
%   Literature:
%     Henkelman RM. Med Phys 1985
%
%   See also SHAPEBRUKERHPMR, RECOBRUKERKSPACE, NSCALE
%
%   08/2019, Keith Michel

%% Parse inputs
if nargin<2,     help(mfilename); return, end
if nargin<3,     exc = []; end
if isempty(exc), exc = 90; end
imgIn = abs(imgIn);
noi   = abs(noi(:));
nStd  = std(noi(:));
[~,~,~,nExc,~] = size(imgIn);
if isscalar(exc), exc = exc*ones(1, nExc); end

%% Check that noise follows Rayleigh distribution
nIcv = mean(noi) / nStd;
if nIcv < 1.8 || nIcv > 2
    warning('calcMzImg:noiseDist', ['Noise might not be true noise. ', ...
        'Inverse CV is %g (should be %g).'], nIcv, sqrt(pi/(4-pi)))
end

%% Rescale to SNR
imgIn = sqrt(2-pi/2) * imgIn ./ nStd;

%% Rescale with excitation angle
[nX,nY,nSlice,nExc,nRep]  = size(imgIn);
imgOut = imgIn ./ repmat(reshape(sind(exc), 1, 1, 1, nExc), ...
    [nX, nY, nSlice, 1, nRep]);
