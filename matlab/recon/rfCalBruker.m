function [p90, bp90, mdl, gf] = rfCalBruker(inDir, saveFiles, noi)
%RFCALBRUKER attempts to calibrate the excitation power for a Bruker hpMR
%   or powcalSINGLEPULSE experiment acquired with variable RF power and
%   either EPI or spectral readouts (for hpMR use an arbitrarily large
%   FOV [>999 cm] to null phase and readout gradients in normal readout mode)
%
%   Usage: [p90, bp90, mdl, gf] = rfCalBruker(inDir, [saveFiles, noi])
%
%   Inputs:
%     inDir     - Integer or char of hpMR scan directory
%     saveFiles - If true save figures and workspace to scan directory
%                   (default = false)
%     noi       - Complex acquisition noise
%                   (only used if not acquired in hpMR scan)
%
%   Outputs:
%     p90  - Power for 90 deg excitation with applied pulse
%     bp90 - Power for 90 deg excitation with 1 ms hard pulse
%     mdl  - Fit object describing model (output from FIT)
%     gf   - Goodness of fit structure (output from FIT)
%
% 09/2019, Keith Michel

funVersion = '20191118';

%% Parse inputs
if ~nargin,             inDir = uigetdir(); end
if isscalar(inDir),     inDir = num2str(inDir); end
if ~exist(inDir, 'dir')
    error('Specified scan directory (%s) does not exist', inDir)
end
if nargin<2,            saveFiles = false; end
if isdeployed
    if nargin<2,        saveFiles   = true; end
end
if nargin<3,            noi       = []; end

if isdeployed
    diary(fullfile(inDir, sprintf('%s-log_%s.txt', mfilename, datestr(now, 30))));
    fprintf('Entering %s at %s\n', mfilename, datestr(now))
    funVersion
end

if strcmp(saveFiles, '0'),  saveFiles = false; end
if ~exist('export_fig', 'file')
    warning(sprintf('%s:export_fig', mfilename), ...
        'export_fig not found, PRINT will be used for exporting figures')
end

%% Read headers
method = readBrukerHeader(fullfile(inDir,'method'));
if ~isempty(regexp(method.Method, 'hpMR', 'once'))
    if strcmpi(method.KAM_VFAOn, 'no')
        error('RF calibration hpMR scan must be run with variable flip angles')
    end
    if strcmpi(method.KAM_FidEcho, 'yes')
        error('RF calibration hpMR scan cannot be run with a fid echo')
    elseif strcmpi(method.KAM_SatEcho, 'yes')
        error('RF calibration hpMR scan cannot be run with a saturation echo')
    end
    if method.KAM_IdealEchoShift > 1e-3
        warning(sprintf('%s:idealEchoShift', mfilename), ...
            'IDEAL echo shift should not be used in RF calibration scan')
    end
    if ~isempty(regexp(lower(method.KAM_ReadMode), 'epi', 'once'))
        fprintf('Performing RF calibration on EPI images\n')
        useIms = true;
    elseif all(method.PVM_FovCm > 999)
        fprintf('Performing RF calibration on spectral readouts\n')
        useIms = false;
    end
    pWatts = method.KAM_VFAWatts;
    if ~strcmpi(method.KAM_ExcMode, 'normal_exc')
        pname = regexprep(method.KAM_ExcMode, '_', ' ');
    else
        pname = method.ExcPulse1Enum;
    end
elseif ~isempty(regexp(method.Method, 'powcalSINGLEPULSE', 'once'))
    useIms = false;
    pWatts = method.MDA_RFPowList;
    pname  = method.ExcPulse1Enum;
else
    error('Not a hpMR or powcalsinglepulse experiment')
end


%% Get signal as function of RF power
[m,k,~,n] = recoBrukerKspace(inDir, [], [], 0);
if useIms
    [~,~,nSlice,nEcho,~] = size(m);
else
    [~,nEcho,nSlice,~]   = size(k);
    nEcho = nEcho - nnz(pWatts < 1e-3);
end
if nSlice > 1
    warning(sprintf('%s:multiSlice', mfilename), ...
        'Using center slice of multislice hpMR scan')
end
if ~isempty(noi)
    fprintf('Using noise provided as input to %s\n', mfilename)
elseif ~isempty(n)
    noi = n;
    fprintf('Using noise echoes acquired in RF calibration scan\n')
elseif useIms == false && any(pWatts < 1e-3)
    noi = k(:,pWatts<1e-3,ceil(nSlice/2),:);
    fprintf('Using noise echoes acquired in RF calibration scan\n')
else
    fprintf('No noise acquired or input, using raw magnitude\n')
end
snrTh = 5;
mSnr  = [];
if useIms
    % Use magnitude of brightest voxel
    m      = sum(m(:,:,ceil(nSlice/2),:,:), 5);
    useSig = true(nEcho, 1);
    if ~isempty(noi)
        n = ifftdim(reshape(noi, size(m, 1), size(m, 2), []), 1:2);
        mSnr    = squeeze(calcMzImg(m, n));
        [~,idx] = max(abs(mSnr(:)));
        [x,y,~] = ind2sub(size(mSnr), idx);
        mSnr    = squeeze(mSnr(x,y,:));
        useSig  = squeeze(mSnr >= snrTh);
        mSnr    = mSnr(useSig);
    end
    m       = squeeze(m);
    [~,idx] = max(abs(m(:)));
    m       = m .* exp(-1i * angle(m(idx)));
    [x,y,~] = ind2sub(size(m), idx);
    sig     = real(squeeze(m(x,y,:)));
else
    % Use fwhm integral of magnitude spectrum
    sz  = size(k);
    sz1 = 2 ^ nextpow2(sz(1));
    k   = [k; zeros([sz1-sz(1), sz(2:end)])];
    m   = ifftdim(sum(k(:,pWatts>=1e-3,ceil(nSlice/2),:), 4), 1);
    useSig = true(nEcho, 1);
    if ~isempty(noi)
        n      = ifftdim(reshape(noi, size(m, 1), []), 1);
        mSnr   = calcMzImg(reshape(m, sqrt(sz1), sqrt(sz1), 1, []), n);
        mSnr   = squeeze(reshape(mSnr, sz1, []));
        mSnr   = max(mSnr, [], 1);
        useSig = mSnr >= snrTh;
        mSnr   = mSnr(useSig);
    end
    m       = squeeze(m);
    [~,idx] = max(abs(m(:)));
    m       = m .* exp(-1i * angle(m(idx)));
    [~,y]   = ind2sub(size(m), idx);
    [~,x]   = hhfw(real(m(:,y)));
    sig     = sum(real(m(ceil(x(1)):floor(x(2)),:)), 1);
end
if ~isempty(noi)
    if ~any(useSig)
        error('All echoes have SNR > %g', snrTh)
    else
        fprintf('Removing %d (of %d) echoes with SNR < %g for RF power calibration\n', ...
            nnz(~useSig), nEcho, snrTh)
    end
end
pWatts = pWatts(pWatts >= 1e-3);
pWatts = pWatts(useSig);
sig    = sig(useSig);

%% Fit the observed real signal as function of RF power
ft        = fittype(@(a, b, x) a * sind(90 * x / b));
xData     = sqrt(50 * pWatts(:));    % W -> V (assuming 50 ohms)
yData     = sig(:);
[mx,idx]  = max(yData);
% if isempty(mSnr)
    [md1,gf1] = fit(xData, yData, ft, 'StartPoint', [mx, xData(idx)]);
% else
%     fprintf('Performing SNR-weighted fit\n')
%     [md1,gf1] = fit(xData, yData, ft, 'StartPoint', [mx, xData(idx)], ...
%         'Weights', mSnr);
% end
% run fit again with inverted signal curve, accept fit with lower sse
yData     = -sig(:);
[mx,idx]  = max(yData);
% if isempty(mSnr)
    [md2,gf2] = fit(xData, yData, ft, 'StartPoint', [mx, xData(idx)]);
% else
%     [md2,gf2] = fit(xData, yData, ft, 'StartPoint', [mx, xData(idx)], ...
%         'Weights', mSnr);
% end
if gf1.sse < gf2.sse
    yData = sig(:);
    mdl   = md1
    gf    = gf1;
else
    mdl   = md2
    gf    = gf2;
end
cvals = coeffvalues(mdl);
p90   = cvals(2)^2 / 50;    % V -> W (assuming 50 ohms)

% use applied pulse length and power integration factor to calculate bp90
pLen = regexp(method.ExcPulse1, '^\((\S+),', 'tokens');
pLen = str2double(cell2mat(pLen{1}));
sInt = regexp(method.ExcPulse1, '^\(\S+,\s+\S+,\s+\S+,\s+\S+,\s+\S+,\s+\S+,\s+(\S*),', 'tokens');
sInt = str2double(cell2mat(sInt{1}));
bp90 = p90 * (sInt * pLen)^2;

fprintf('R^2 = %.06f\n', gf.rsquare)
fprintf('p90 = %.04f W (Applied pulse)\nref90 = %.04f W (1 ms Hard pulse)\n', ...
    p90, bp90)

%% Plot results
exca = 90 * xData / cvals(2);   % calculate actual acquired excitation angles

if ~isempty(regexp(method.Method, 'hpMR', 'once')) && ...
        strcmpi(method.PVM_DeriveGains, 'yes')
    fig = figure('Color', 'w', 'Position', 100*[1 1 12 5], ...
        'Name', sprintf('%s - p90=%.04fW', mfilename, p90));
    subplot 122
    pDeg = method.KAM_VFADegrees(method.KAM_VFAWatts >= 1e-3);
    pDeg = pDeg(useSig);
    plot(pDeg, exca, 'o', 'MarkerSize', 9, 'MarkerEdgeColor', ...
        [0.2 0.6275 0.1725], 'MarkerFaceColor', [0.698 0.8745 0.5412]);
    hold on
    plot([0 90], [0 90], 'LineWidth', 1.5, 'Color', [0.2 0.6275 0.1725]);
    xlabel('Nominal Excitation Angle (deg)')
    ylabel('Measured Excitation Angle (deg)')
    set(gca, 'fontsize', 14)
    grid on
    subplot 121
else
    fig = figure('Color', 'w', 'Name', sprintf('%s - p90=%.04fW', mfilename, p90));
end

plot(pWatts, yData, 'o', 'MarkerSize', 9, 'MarkerEdgeColor', ...
    [0.1216 0.4706 0.7059], 'MarkerFaceColor', [0.6510 0.8078 0.8902]);
hold on
excs = linspace(min(exca), max(exca), 512);
pows = (excs * cvals(2) / 90).^2 / 50;
plot(pows, cvals(1)*sind(excs), 'LineWidth', 1.5, ...
    'Color', [0.1216 0.4706 0.7059])
xlabel('RF Power (W)')
if useIms
    ylabel('Real Image Signal')
else
    ylabel('Real Signal')
end
set(gca, 'fontsize', 14)
title(pname, 'fontsize', 15);
grid on

text(0.05, 0.85, sprintf('R^2 = %.04f', gf.rsquare), 'Units', 'normalized', ...
    'Fontsize', 14);
text(0.6, 0.15, {sprintf('p90 = %.04f W (Applied pulse)', p90), ...
    sprintf('ref90 = %.04f W (1 ms Hard pulse)', bp90)}, ...
    'Units', 'normalized', 'Fontsize', 14, 'HorizontalAlignment', 'center');

if saveFiles
    pngFile = fullfile(inDir, sprintf('%s-%s', mfilename, datestr(now, 30)));
    if exist('export_fig', 'file')
        export_fig(fig, pngFile, '-a1')
    else
        print(fig, pngFile, '-dpng')
    end
    % save(fullfile(inDir, sprintf('%s-%s.mat', mfilename, datestr(now, 30))), ...
    %     'inDir', 'noi', 'pname', 'p90', 'bp90', 'sInt', 'pLen', 'method', ...
    %     'mdl', 'gf', 'mSnr', 'funVersion', 'pWatts', 'sig', 'use*');
    saveVars = {'inDir', 'noi', 'pname', 'p90', 'bp90', 'sInt', 'pLen', ...
        'method', 'mdl', 'gf', 'mSnr', 'funVersion', 'pWatts', 'sig', ...
        'useIms', 'useSig'};
    uisave(saveVars, fullfile(inDir, ...
        sprintf('%s-workspace_%s.mat', mfilename, datestr(now, 30))));
end

if isdeployed
    waitfor(gcf);
    fprintf('Leaving %s at %s\n', mfilename, datestr(now))
end

end

