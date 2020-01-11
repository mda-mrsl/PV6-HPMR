function ksp = ghostCorrEpi(ksp_in)
%GHOSTCORREPI applies heuristic image-based N/2 ghost correction for EPI 
%
%   Usage: ksp = ghostCorrEpi(ksp_in)
%
%       where ksp_in is complex EPI k-space data reshaped to size 
%               [nRead, nPhase, nSlice, nEcho, nRep] as by shapeBrukerHpmr
%             ksp is ghost-corrected EPI k-space data 
%
%   Literature: 
%     Buonocore MH, Gao L. MRM 1997
% 
%   See also SHAPEBRUKERHPMR, RECOBRUKERKSPACE, CREATE2DEPI
%
%   09/2019, Keith Michel

funVersion = 'v20191025';

%% Check inputs
if nargin<1,	help(mfilename); return; end
[mtx,mty,nSlice,nEcho,nRep] = size(ksp_in);
if mtx~=mty
    error('ghostCorrEpi:mtx', 'Non-square matrices have not been verified')
end

if isdeployed
    fprintf(1, 'Entering %s (%s) at %s\n', mfilename, funVersion, datestr(now));
end

%% Heuristic image-based N/2 ghost correction for each echo
np   = 40;                  % number of 0th and 1st order phase coeffs to test
p0mx = pi/2;                % +/- max for 0th order phase coeff
p1mx = pi/8;                % +/- max for 1st order phase coeff
p0   = linspace(-p0mx, p0mx, np);
p1   = linspace(-p1mx, p1mx, np);
dist = floor(sqrt(mty))-1;  % size of ghost region at edge of image matrix
mask = zeros(mtx, mty);
mask(1:mtx,1:dist) = 1;
mask = mask + fliplr(mask);
uSl  = ceil(nSlice/2);      % use center slice
gMag = zeros(np, np, nEcho);
gcIm = zeros(size(ksp_in));
for ii=1:nEcho
    kOdd = zeros(mtx, mty);
    kEvn = zeros(mtx, mty);
    
    % use repetition with maximum signal
    curKsp     = squeeze(ksp_in(:,:,uSl,ii,:));
    [~,idx]    = max(abs(curKsp(:)));
    [~,~,uRep] = ind2sub(size(curKsp), idx);
    kOdd(:,1:2:end) = curKsp(:,1:2:end,uRep);
    kEvn(:,2:2:end) = flipud(curKsp(:,2:2:end,uRep));
    mOdd = ifftdim(kOdd);
    mEvn = ifftdim(kEvn);
    
    % locate phase coeffs that minimize signal magnitude in ghost regions
    for jj=1:np^2
        [a,b] = ind2sub([np,np], jj);
        pmat  = repmat(p0(a) * linspace(-p1(b), p1(b), mtx).', [1,mty]);
        tmp   = abs(mOdd .* exp(1i*pmat) + mEvn .* exp(-1i*pmat));
        gMag(a,b,ii) = sum(sum(tmp .* mask));
    end
    [~,idx] = min(gMag(a,b,ii));
    [a,b]   = ind2sub([np,np], idx);
    pmat    = repmat(p0(a) * linspace(-p1(b), p1(b), mtx).', ...
        [1,mty,nSlice,1,nRep]);
    
    % Recon images using optimal phase coefficients
    kOdd = zeros(mtx, mty, nSlice, 1, nRep);
    kEvn = zeros(mtx, mty, nSlice, 1, nRep);
    kOdd(:,1:2:end,:,:,:) = ksp_in(:,1:2:end,:,ii,:);
    kEvn(:,2:2:end,:,:,:) = flipud(ksp_in(:,2:2:end,:,ii,:));
    mOdd = ifftdim(kOdd, 1:2);
    mEvn = ifftdim(kEvn, 1:2);
    gcIm(:,:,:,ii,:) = mOdd.*exp(1i*pmat) + mEvn.*exp(-1i*pmat);
end

% Return k-space of ghost-corrected images
ksp = fftdim(gcIm, 1:2);

if isdeployed
    fprintf(1, 'Leaving %s at %s\n', mfilename, datestr(now));
end