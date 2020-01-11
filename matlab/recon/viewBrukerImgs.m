function [imgs, reco, acqp] = viewBrukerImgs(dir, recoNumber, showImgs)
% VIEWBRUKERIMGS reads in Paravision reconstructed images and displays as montage
%       For spectroscopy, the spectra are plotted
%
%   Usage: [imgs, reco, acqp] = viewBrukerImgs(dir, [recoNumber, showImgs])
%
%       where dir is a Paravision scan directory
%               if omitted, working directory is used
%             recoNumber is the PV reconstruction number
%               if omitted, 1 [i.e. imgs in fullfile(dir, 'pdata/ *1* /2dseq')]
%             showImgs indicates whether to display output as montage
%               if omitted, true
%             imgs is the array of Paravision reconstructed images
%             reco is the Paravision reco header file structure
%             acqp is the Paravision acqp header file structure
%
%   See also IMAGESC, SHOWMONTAGE, READBRUKERSTUDY, RECOBRUKERKSPACE
%
%   06/2019, Keith Michel

funVersion = 'v20191007';

%%
if ~nargin
    if isdeployed, dir = uigetdir('', 'Select Bruker Experiment Directory');
    else,          dir = '.'; end
end
dir = num2str(dir);
if nargin<2, recoNumber = 1; end
if nargin<3, showImgs = true; end
recoNumber = num2str(recoNumber);
fullDir    = fullfile(dir, 'pdata', recoNumber);

if isdeployed
    if ~regexpi(mfilename, 'viewbrukerimages', 'once')
        diary(fullfile(dir, sprintf('%s-log_%s.txt', mfilename, datestr(now, 30))));
    end
    fprintf(1, 'Entering %s (%s) at %s\n', mfilename, funVersion, datestr(now));
end

%%
acqp = readBrukerHeader(fullfile(dir, 'acqp'));
reco = readBrukerHeader(fullfile(fullDir, 'reco'));
f    = fopen(fullfile(fullDir, '2dseq'));
if strcmpi(reco.RECO_wordtype, '_16BIT_SGN_INT')
    data = fread(f, inf, 'int16');
elseif strcmpi(reco.RECO_wordtype, '_32BIT_SGN_INT')
    data = fread(f, inf, 'int32');
else
    error('viewBrukerImgs:RECO_wordtype', 'Unrecognized image data type')
end
fclose(f);

%%
imgSize = reco.RECO_size;
offset  = reco.RECO_map_offset;
slope   = reco.RECO_map_slope;
if isfield(reco, 'RecoNumRepetitions') % PV 6
    nRep = reco.RecoNumRepetitions;
    nImg = reco.RecoObjectsPerRepetition;
else
    nRep = acqp.NR;
    nImg = numel(data) / prod(imgSize) / nRep;
end
if numel(imgSize) > 1
    imgSize = imgSize(1:2);
    rot  = reco.RECO_transposition;
    nXY  = max(imgSize);
    nZ   = numel(data)/prod(imgSize);
    imgs = nan(nXY, nXY, nZ);
    
    data = reshape(data, imgSize(1), imgSize(2), []);
    if numel(rot) < nZ
        rot    = repmat(rot, 1, nZ/numel(rot));
    end
    if numel(offset) < nZ
        offset = repmat(offset, 1, nZ/numel(offset));
        slope  = repmat(slope, 1, nZ/numel(slope));
    end
    for ii = 1:nZ
        if rot(ii)
            imgs(1:imgSize(2),1:imgSize(1),ii) = ...
                fliplr(reshape((data(:,:,ii) + offset(ii)) / slope(ii), ...
                imgSize(2), imgSize(1)));
        else
            imgs(1:imgSize(2),1:imgSize(1),ii) = ...
                rot90(fliplr((data(:,:,ii) + offset(ii)) / slope(ii)), 1);
        end
    end
    try
        imgs = reshape(imgs, nXY, nXY, nImg, nRep);
    catch
        warning('viewBrukerImgs:imgReshape', ...
            '4D reshape failed, returning 3D image array')
        imgs = reshape(imgs, nXY, nXY, []);
    end
    uRot  = unique(rot);
    uSize = unique(imgSize);
    if numel(uSize) > 1 && numel(uRot) == 1
        if imgSize(1) > imgSize(2)
            imgs(imgSize(2)+1:end,:,:,:) = [];
        else
            imgs(:,imgSize(1)+1:end,:,:) = [];
        end
    end
    if showImgs
        showMontage(imgs);
    end
else
    imgs = (data + offset) / slope;
    if numel(imgs) > imgSize
        imgs = reshape(imgs, imgSize, []);
    end
    if showImgs
        plotSpec(imgs, reco.RECO_sw);
    end
end
if showImgs
    set(gcf, 'Name', sprintf('%s Images - %s', mfilename, acqp.ACQ_scan_name))
end

if isdeployed
    if ~regexpi(mfilename, 'viewbrukerimages', 'once')
        % save(fullfile(dir, sprintf('%s-workspace_%s', ...
        %     mfilename, datestr(now, 30))), 'imgs', 'acqp', 'reco')
        saveVars = {'imgs', 'acqp', 'reco'};
        uisave(saveVars, fullfile(dir, sprintf('%s-workspace_%s', ...
            mfilename, datestr(now, 30))));
        if ~ismatrix(imgs) && showImgs
            viewover(squeeze(img), acqp.ACQ_scan_name)
            waitfor(gcf)
        end
    end
    fprintf(1, 'Leaving %s at %s\n', mfilename, datestr(now));
end

