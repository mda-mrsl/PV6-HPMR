function [ims, bkg, noi] = hpSnapEpi(inDir, saveFiles, fuseIms, bkgDir, noi, idealHz)
%HPSNAPEPI reconstructs images from hyperpolarized [1-13C]pyruvate hpMR
%        experiment acquired with EPI readouts and either IDEAL encoding or
%        spectral-spatial single-metabolite excitations.
%
%   Usage: [ims, bkg, noi] = hpSnapEpi(inDir, [saveFiles, fuseIms, bkgDir, noi])
%
%   Inputs:
%     inDir     - Integer or char of hpMR scan directory
%     saveFiles - If true save figures and workspace to scan directory
%                 If 2 and fuseIms is false, hpFuseIms is not called
%                   (default = false, 2 ifdeployed)
%     fuseIms   - If true call hpFuseIms to overlay images
%                   (default = false)
%     bkgDir    - Integer or char of background scan directory
%     noi       - Complex acquisition noise, *not* FT'ed
%                   (if not provided, noise acquired in hpMR scan is used)
%
%   Outputs:
%     ims - 5D array of dynamic magnitude images
%           size: [mtx, mty, nSlices, nEchoes, nReps]
%     bkg - Background images (directly from hpFuseIms)
%     noi - Complex acquisition noise used for SNR calculation, if present
%
% 08/2019, Keith Michel

funVersion = 'v20191106';

%% Parse inputs
if ~nargin,             inDir = uigetdir(); end
if isscalar(inDir),     inDir = num2str(inDir); end
if ~exist(inDir, 'dir')
    error('Specified scan directory (%s) does not exist', inDir)
end
if nargin<2,            saveFiles = false; end
if isdeployed
    if nargin<2,        saveFiles = 2; end
end
if nargin<3,            fuseIms   = false; end
if isdeployed
    if nargin<3,        fuseIms   = true; end
end
if nargin<4,            bkgDir    = []; end
if nargin<5,            noi       = []; end
if nargin<6,            idealHz   = []; end

if isdeployed
    diary(fullfile(inDir, sprintf('%s-log_%s.txt', mfilename, datestr(now, 30))));
    fprintf(1, 'Entering %s (%s) at %s\n', mfilename, funVersion, datestr(now));
end

if strcmp(saveFiles, '0'),  saveFiles = false; end
if strcmp(fuseIms, '0'),  	fuseIms = false; end
if ~exist('export_fig', 'file')
    warning('hpSnapEpi:export_fig', ...
        'export_fig not found, PRINT and GIF will be used for exporting figures')
end

%% Read headers
method = readBrukerHeader(fullfile(inDir, 'method'));
if isempty(regexp(method.Method, 'hpMR', 'once'))
    error('Not a hpMR experiment')
end
if strcmpi(method.KAM_FidEcho, 'yes')
    params = struct();
    if isdeployed && verLessThan('matlab', '8.4'), params.showPlots = 2;
    else,                                         params.showPlots = true; end
    if ~saveFiles && ~fuseIms,                    params.showPlots = false; end
    hpSliceSpectra = hpFidEcho(inDir, params);
end
if isempty(regexp(lower(method.KAM_ReadMode), 'epi', 'once')) && ...
    isempty(regexp(lower(method.KAM_ReadMode), 'flyback', 'once'))
    error('Specified scan does not use an EPI readout.');
end
if method.KAM_IdealEchoShift > 1e-3
    idealOn = true;
else%if  ~isempty(regexp(lower(method.KAM_ExcMode), 'spsp', 'once'))
    idealOn = false;
end

%% hpMR recon
[m,k,~,n,~,s] = recoBrukerKspace(inDir, [], [], 2, idealHz);
[mtx,mty,nSlice,nEcho,nRep] = size(m);
sImg = [];
if ~isempty(s)
    sImg  = ifftdim(reshape(s, mtx, mty, []), 1:2);
    sFreq = [];
    if isfield(method, 'KAM_SatFreq')
        sFreq = method.KAM_SatFreq;
    end
end
if ~isempty(noi)
    fprintf(1, 'Using noise provided as input to %s\n', mfilename);
elseif ~isempty(n)
    noi  = n;
    fprintf(1, 'Using noise echoes acquired in hpMR scan\n');
end
if numel(noi) >= mtx*mty
    noi  = noi(1:mtx*mty*floor(numel(noi)/mtx/mty));
    nn   = reshape(noi, mtx*mty, floor(numel(noi)/mtx/mty));
    nn   = ifftdim(nn, 1:2);
    ims  = calcMzImg(m, nn);
    sImg = calcMzImg(sImg, nn);
    fprintf(1, 'Returning SNR images (noise value = %g, inverse CV = %g)\n', ...
        std(abs(nn(:)))/sqrt(2-pi/2), mean(abs(nn(:)))/std(abs(nn(:))));
    useSnr = true;
elseif ~isempty(noi)
    warning('hpSnapEpi:inputNoiseSize', ['Noise array provided as input ', ...
        'to hpSnapEpi is not large enough to apply inverse FT.'])
    ims  = abs(m);
    sImg = abs(sImg);
    fprintf(1, 'Bad noise input, returning raw magnitude images\n');
    useSnr = false;
else
    ims  = abs(m);
    sImg = abs(sImg);
    fprintf(1, 'No noise acquired or input, returning raw magnitude images\n');
    useSnr = false;
end
clear m n s

tRange  = 5:55;  % Repetitions to average in pngs and show in gifs for hpFuseIms
if ~idealOn
    if strcmpi(method.KAM_SpecExc, 'yes')
        txFreqs = method.KAM_VOFTxHertz;
    else
        txFreqs = zeros(1, nEcho);
    end

    vfaDeg = method.KAM_VFADegrees;
    vfaDeg = vfaDeg(vfaDeg > 0.1);
    
    pyrImgs = squeeze(ims(:,:,:,2,:));
    lacImgs = squeeze(ims(:,:,:,1,:));
    ureImgs = sImg;
else
    pyrImgs = squeeze(ims(:,:,:,1,:));
    lacImgs = squeeze(ims(:,:,:,2,:));
    ureImgs = squeeze(ims(:,:,:,4,:));
    txFreqs = {'Pyruvate', 'Lactate', 'Hydrate', 'Urea', 'Alanine', 'Bicarbonate'};
end

savePngs = false;
saveGifs = false;
if saveFiles == 1
    savePngs = true;
    saveGifs = true;
end
if saveFiles == 1 || fuseIms
    save(fullfile(inDir, sprintf('%s-workspace', mfilename)))
elseif isdeployed
    saveVars = {'bkgDir', 'hpSliceSpectra', 'funVersion', 'fuseIms', ...
        'idealOn', 'ims', 'inDir', 'k', 'lacImgs', 'method', 'mtx', 'mty', ...
        'nEcho', 'nRep', 'nSlice', 'noi', 'pyrImgs', 'sImg', 'saveFiles', ...
        'saveGifs', 'savePngs', 'tRange', 'txFreqs', 'ureImgs', 'useSnr'};
    uisave(saveVars, ...
        fullfile(inDir, sprintf('%s-workspace', mfilename)));
end

%% Show images
if ~isempty(sImg)
    showMontage(sImg);
    set(gcf, 'Color', 'w', 'Name', 'hpMR Images - Saturation Echo');
    title(sprintf('Saturation Echo Images, %+g Hz Tx Offset', sFreq), 'FontSize', 15);
    h = colorbar;
    set(h, 'FontSize', 10);
    if useSnr
        set(get(h, 'Title'), 'string', 'SNR')
    end
    if saveFiles > 0
        pngFile = fullfile(inDir, sprintf( ...
            '%s-dynamicMontage-SaturationEcho', mfilename));
        if exist('export_fig', 'file')
            export_fig(gcf, pngFile, '-a1')
        else
            print(gcf, pngFile, '-dpng')
        end
    end
end
for ii = 1:nEcho
    if ~idealOn
        txt = sprintf('%+g Hz Tx Offset', txFreqs(ii));
        if strcmp(method.KAM_VFAOn, 'Yes')
            txt = [txt, sprintf(', %g Degrees', vfaDeg(ii))];
        end
    else
        txt = ['IDEAL ', txFreqs{ii}];
    end
    for jj = 1:nSlice
        showMontage(ims(:,:,jj,ii,:));
        set(gcf, 'Color', 'w', 'Name', ['hpMR Images - ', txt]);
        title(txt, 'FontSize', 15);
        h = colorbar;
        set(h, 'FontSize', 10);
        if useSnr
            set(get(h, 'Title'), 'string', 'SNR')
        end
        if saveFiles > 0
            pngFile = fullfile(inDir, sprintf( ...
                    '%s-dynamicMontage-Echo%d-Slice%d', mfilename, ii, jj));
            if exist('export_fig', 'file')
                export_fig(gcf, pngFile, '-a1')
            else
                print(gcf, pngFile, '-dpng')
            end
        end
    end
end

if saveFiles > 0
    % Export residual spectra and histogram figures
    figs  = findobj('Type', 'figure');
    for ii = 1:numel(figs)
        pngFile = [];
        switch get(figs(ii), 'Name')
            case 'hpMR Residual Spectra'
                pngFile = fullfile(inDir, sprintf( ...
                    '%s-residualSpectra', mfilename));
            case 'hpMR Noise'
                pngFile = fullfile(inDir, sprintf( ...
                    '%s-noiseHistogram', mfilename));
        end
        if ~isempty(pngFile)
            if exist('export_fig', 'file')
                export_fig(figs(ii), pngFile, '-a1')
            else
                print(figs(ii), pngFile, '-dpng')
            end
        end
    end
end

bkg = [];
if fuseIms
    bkg = hpFuseIms(inDir, bkgDir, savePngs, saveGifs, true);
elseif saveFiles == 1
    bkg = hpFuseIms(inDir, bkgDir, savePngs, saveGifs, false);
end

if isdeployed
    if fuseIms
        waitfor(gcf);
    else
        close all
    end
    fprintf(1, 'Leaving %s at %s\n', mfilename, datestr(now));
end

