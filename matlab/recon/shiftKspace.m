function ksp = shiftKspace(ksp_in, shift)
%SHIFTKSPACE applies phase shift to k-space data
%
%   Usage: ksp = shiftKspace(ksp_in, shift)
%
%       where ksp_in is raw data k-space data
%             shift is the desired shift as fraction of image domain FOV 
%               can be a vector with shift value for each dimension of ksp
%
%   See also RECOBRUKERKSPACE, SHAPEBRUKERHPMR
%
%   06/2019, Keith Michel

%% Parse inputs
if ~nargin,     help(mfilename); return; end
if numel(shift) > ndims(ksp_in)
    error('shiftKspace:shiftDim', ...
        'Shift requested in %d dimensions but ksp_in is %d-D', ...
            numel(shift), ndims(ksp_in))
end

%% Apply linear phase shift to k-space data
sz    = size(ksp_in);
shift = [shift(:)', zeros(1, numel(sz)-numel(shift))];
ksp   = ksp_in;
for ii = 1:numel(sz)
    p   = 1i*2*pi*shift(ii)*(1:sz(1));
    ksp = ksp .* exp(repmat(p.', [1, sz(2:end)]));
    ksp = shiftdim(ksp, 1);
    sz  = size(ksp);
end