function bkg = hpFuseIms(hpDir, inDir, savePngs, saveGifs, showOver, saveDcms)
%HPFUSEIMS overlays images from hpMR experiment onto FLASH or RARE 1H background images
%
%   Usage: bkg = hpFuseIms(hpDir, inDir, [saveJpgs, saveGifs, showOver, saveDcms])
%
%   Inputs:
%     hpDir     - Integer or char of hpMR scan directory
%     inDir     - Integer or char of FLASH or RARE scan directory
%     savePngs  - If true save pngs to hpMR directory
%                   (default = false)
%     saveGifs  - If true save gifs to hpMR directory
%                   (default = false)
%     showOver  - If true show viewover figure with fused imgs
%                   (default = true)
%     saveDcms  - If true save DICOMs to study directory
%                   (default = true)
%
%   Output:
%     bkg - Background images
%
%
% 08/2019, Keith Michel

funVersion = 'v20191217';
if isdeployed
    fprintf(1, 'Entering %s (%s) at %s\n', mfilename, funVersion, datestr(now));
end

%% Parse inputs
if ~nargin,             hpDir = []; end
if isempty(hpDir),      hpDir = uigetdir('', 'Select hpMR folder');
elseif isscalar(hpDir), hpDir = num2str(hpDir); end
if nargin<2,            inDir = []; end
if isempty(inDir),      inDir = uigetdir(fileparts(hpDir), 'Select underlay folder');
elseif isscalar(inDir), inDir = num2str(inDir); end
if nargin<3,            savePngs = false; end
if nargin<4,            saveGifs = false; end
if nargin<5,            showOver = true; end
if nargin<6,            saveDcms = true; end

%% Read in hpMR data and method header
hpFile = fullfile(hpDir, 'hpSnapEpi-workspace.mat');
if exist(hpFile, 'file')
    fprintf(1, 'Reading hpMR image data for scan %s\n', fileparts(hpDir));
    load(hpFile, 'pyrImgs', 'lacImgs', 'ureImgs', 'method', 'nSlice', ...
        'nEcho', 'nRep', 'useSnr', 'idealOn', 'tRange')
else
    error('hpMR experiment must be reconstructed with hpSnapEpi first')
end
hpMethod = method;
if isempty(ureImgs)
    ureImgs = zeros(size(pyrImgs));
end
for ii = 1:nRep  % The following is backwards compatible with R2013a
    pyrImgs(:,:,ii) = flipud(rot90(pyrImgs(:,:,ii)));
    lacImgs(:,:,ii) = flipud(rot90(lacImgs(:,:,ii)));
    ureImgs(:,:,ii) = flipud(rot90(ureImgs(:,:,ii)));
end

%% Read in background image and method header
method = readBrukerHeader(fullfile(inDir,'method'));
if ~(strcmpi(method.Method, 'bruker:rare') || strcmpi(method.Method, 'bruker:flash'))
    error('Background image not a FLASH or RARE experiment')
end
% TODO: additional checks for consistent scan geometries

bkg = recoBrukerKspace(inDir, [], [], false);
bkg = flipdim(bkg, 3);  % spec-spat slices proceed from head to foot
if size(bkg, 3) == 8    % Hard code for Aug 2019 Lai expts
    bkg = nscale(sum(abs(bkg(:,:,3:6)), 3));
else
    bkg = nscale(abs(bkg(:,:,ceil(size(bkg,3)/2)))); % Center slice
end
bkg = flipud(rot90(bkg));

%% Export DICOMs
if saveDcms
    subject = readBrukerHeader(fullfile(fileparts(hpDir), 'subject'));

    acqp               = readBrukerHeader(fullfile(inDir, 'acqp'));
    reco               = readBrukerHeader(fullfile(inDir, 'pdata', '1', 'reco'));
    acqp.ACQ_scan_name = [acqp.ACQ_scan_name, ' - hpBackground'];
    exportImgs(bkg, 'dicom', fileparts(hpDir), subject, method, acqp, reco);

    hpAcqp  = readBrukerHeader(fullfile(hpDir, 'acqp'));
    hpReco  = readBrukerHeader(fullfile(hpDir, 'pdata', '1', 'reco'));
    % TODO: correct hpAcqp.ACQ_flip_angle for pyr, lac, urea
    hpMethod.PVM_SpatResol = hpMethod.PVM_Fov ./ size(pyrImgs, 1);
    hpScanName             = hpAcqp.ACQ_scan_name;
    hpAcqp.ACQ_scan_name   = [hpScanName, ' - hpPyruvate'];
    exportImgs(pyrImgs, 'dicom', fileparts(hpDir), subject, hpMethod, hpAcqp, hpReco);
    hpAcqp.ACQ_scan_name   = [hpScanName, ' - hpLactate'];
    exportImgs(lacImgs, 'dicom', fileparts(hpDir), subject, hpMethod, hpAcqp, hpReco);
    if nnz(ureImgs)
        hpAcqp.ACQ_scan_name   = [hpScanName, ' - hpUrea'];
        exportImgs(ureImgs, 'dicom', fileparts(hpDir), subject, hpMethod, hpAcqp, hpReco);
    end
end

%% Resize and overlay 13C images
cmap1   = brighten(brewermap(256, '*greens'), -0.7);
cmap2   = brighten(brewermap(256, '*oranges'), -0.7);
cmap3   = brighten(brewermap(256, '*blues'), -0.7);
% tRange  = 5:55;  % Repetitions to average in pngs and show in gifs
if savePngs || saveGifs
    outDir = fullfile(hpDir, mfilename);
    if ~exist(outDir, 'dir')
        mkdir(outDir)
    end
    % The following is backwards compatible with R2013a
    tmpp = zeros(size(bkg,1), size(bkg, 2), nRep);
    tmpl = zeros(size(bkg,1), size(bkg, 2), nRep);
    tmpu = zeros(size(bkg,1), size(bkg, 2), nRep);
    for ii = 1:nRep
        tmpp(:,:,ii) = imresize(pyrImgs(:,:,ii), ...
            flipdim(method.PVM_Matrix, 2), 'bilinear');
        tmpl(:,:,ii) = imresize(lacImgs(:,:,ii), ...
            flipdim(method.PVM_Matrix, 2), 'bilinear');
        tmpu(:,:,ii) = imresize(ureImgs(:,:,ii), ...
            flipdim(method.PVM_Matrix, 2), 'bilinear');
    end
    pyrImgsRsz = tmpp; clear tmpp
    lacImgsRsz = tmpl; clear tmpl
    ureImgsRsz = tmpu; clear tmpu
    % save(fullfile(outDir, sprintf('%s-workspace', mfilename)))
    saveVars = {'bkg', 'hpDir', 'hpMethod', 'funVersion', 'idealOn', ...
        'lacImgs', 'lacImgsRsz', 'method', 'pyrImgs', 'pyrImgsRsz', ...
        'tRange', 'ureImgs', 'ureImgsRsz', 'useSnr'};
    uisave(saveVars, ...
        fullfile(outDir, sprintf('%s-workspace', mfilename)))
else
    if showOver && nnz(ureImgs)
        figure('position', [50 50 800 800], 'Name', 'HP 13C Overlay', ...
            'color', 'w');
        ax1 = axes('position', [0 0 0.98 1]);
        over = mean(ureImgs(:,:,tRange), 3);
        overlayImages(bkg, over, gray(256), cmap1, ax1);
        h2  = colorbar;
        set(h2, 'FontSize', 12)
        if useSnr
            set(get(h2, 'Title'), 'string', 'Mean SNR')
        end
        title(sprintf('Urea, Reps %02d to %02d Averaged', ...
            tRange(1), tRange(end)), 'fontsize', 16)
        viewover([pyrImgs, lacImgs*5], [bkg, bkg], 'title', ...
            'Pyruvate and Lactate(x5)', 'zoom', 4)
        waitfor(gcf)
    else
        figure('position', [50 50 800 800], 'Name', 'HP 13C Overlay', ...
            'color', 'w');
        ax1 = axes('position', [0 0 0.98 1]);
    end
end

if savePngs || saveGifs
    fprintf(1, 'Exporting image overlays...');
    tic
    fig = figure('position', [50 50 800 800], 'Name', 'HP 13C Overlay', ...
        'color', 'w');
    pyrMax = max(pyrImgsRsz(:));
    lacMax = max(lacImgsRsz(:));
    ureMax = max(ureImgsRsz(:));
    mets = {'Urea', 'Lactate', 'Pyruvate'};
    mshs = {'ure', 'lac', 'pyr'};
%     iterRange = 1:3;
%     if nnz(ureImgs)
%         iterRange = 2:3;
%     end
    for ii = 1:3
        met = mets{ii};
        msh = mshs{ii};
        if ~exist(sprintf('%sImgs', msh), 'var') || ~eval(sprintf('nnz(%sImgs)', msh))
            fprintf('Skipping %s Images', met);
            continue
        end
        eval(sprintf('cmap = cmap%d;', ii))

        clf(fig)
        figure(fig)
        ax1 = axes('position', [0 0 0.98 1]);
        eval(sprintf('over = mean(%sImgsRsz(:,:,tRange), 3);', msh))
        overlayImages(bkg, over, gray(256), cmap, ax1);
        h2  = colorbar;
        set(h2, 'FontSize', 12)
        if useSnr
            set(get(h2, 'Title'), 'string', 'Mean SNR')
        end
        title(sprintf('%s, Reps %02d to %02d Averaged', ...
            met, tRange(1), tRange(end)), 'fontsize', 16)
        drawnow
        if savePngs
            pngFile = fullfile(outDir, sprintf('AvgIm-%s', met));
            if exist('export_fig', 'file')
                export_fig(fig, pngFile, '-a1')
            else
                print(fig, pngFile, '-dpng')
            end
        end

        if showOver && ii == 1
            viewover(cat(2, pyrImgs, lacImgs*5), [bkg, bkg], 'title', ...
                'Pyruvate and Lactate(x5)', 'zoom', 4)
            waitfor(gcf)
        end

        if ~strcmp(mshs{ii}, 'ure')
            for jj = tRange
                clf(fig)
                figure(fig)
                ax1 = axes('position', [0 0 0.98 1]);
                eval(sprintf('over = %sImgsRsz(:,:,jj);', msh))
                [~,imh] = overlayImages(bkg, over, gray(256), cmap, ax1);
                eval(sprintf('caxis([0 %sMax]);', msh))
                eval(sprintf('set(imh(2), ''AlphaData'', over ./ %sMax)', msh))
                h2  = colorbar;
                set(h2, 'FontSize', 12)
                if useSnr
                    set(get(h2, 'Title'), 'string', 'SNR')
                end
                title(sprintf('%s, Rep %02d', met, jj), 'fontsize', 16)
                drawnow
                if savePngs
                    pngFile = fullfile(outDir, sprintf( ...
                        'DynIms-%s-Rep%02d', met, jj));
                    if exist('export_fig', 'file')
                        export_fig(fig, pngFile, '-a1')
                    else
                        print(fig, pngFile, '-dpng')
                    end
                end
                if saveGifs
                    gifFile = fullfile(outDir, sprintf('DynGif-%s', met));
                    if jj==tRange(1)
                        gif([gifFile, '.gif'], 'DelayTime', 1/5, ...
                            'frame', fig)
                    elseif jj <= tRange(end)
                        gif
                    end
                    if jj == tRange(end)
                        gif('clear')
                    end
                end
            end
        end
    end
    close(fig);
    fprintf(1, ' Done!\n');
    toc
end

if isdeployed
    fprintf(1, 'Leaving %s at %s\n', mfilename, datestr(now));
end

