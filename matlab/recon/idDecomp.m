function dataOut = idDecomp(dataIn,chemShifts,echoTimes,dwMap,T2s,showPlot)
%IDDECOMP decomposes chemical shift data from complex X- or K-space data at
%         given echo times. Optionally use dwell time correction of K-space
%         data, correct for T2* decay and show sim plots.
%
%   Usage: dataOut = idDecomp(dataIn,chemShifts,echoTimes,[dwMap],[T2s])
%
%   Inputs:
%     dataIn      - Array of complex multiecho data
%     chemShifts  - Vector of chemical shifts (Hz)
%     echoTimes   - Vector of echo times (ms)
%     dwMap       - Vector of pointwise relative acquisition times
%                   in K-space (ms)
%         [OPTIONAL], default  = []
%     T2s         - T2* relaxation time constant (ms)
%         [OPTIONAL], default = inf
%     showPlot    - show IDEAL performance simulation plots (see IDSIMPLOT)
%         [OPTIONAL], default = false
%
%   Output:
%     dataOut     - Array of data for individual chemical shifts
%
%   Literature:
%     Brodsky EK, et al. MRM 2008. DOI 10.1002/mrm.21580
%
% 07/2019, Keith Michel

if nargin<3, help(mfilename); return; end

%% Parse inputs
if ~exist('dwMap','var'),       dwMap = []; end
if ~exist('T2s','var'),         T2s = []; end
if isempty(T2s),                T2s = inf; end
if ~exist('showPlot','var'),    showPlot = []; end
if isempty(showPlot),           showPlot = false; end

dataSize = size(dataIn);
nEchoes = numel(echoTimes);
nShifts = numel(chemShifts);
if numel(dataSize) ~= 2 || dataSize(2) ~= nEchoes
    error('dataIn must be 2D array of size [nPoints nEchoes]')
end

if ~isempty(dwMap)
    dwMap = dwMap(:);
    if numel(dwMap) ~= dataSize(1)
        error('Dwell time map size must match dataIn')
    end
else
    dwMap = zeros(dataSize(1),1);
end

chemShifts = chemShifts(:).';
echoTimes = echoTimes(:);

%% NSA and condition number simulations
i2p = 1i * 2 * pi;
if showPlot
    idSimPlot(chemShifts, diff(echoTimes(1:2)), nEchoes);
end

%% Pointwise IDEAL
phasorMatrix = exp( i2p * 1e-3 * echoTimes * chemShifts ) .* ...
    repmat( exp( -echoTimes / T2s), 1, nShifts);
dataOut = (dataIn * pinv(phasorMatrix).') .* ...
    exp( -i2p * 1e-3 * dwMap * chemShifts);

