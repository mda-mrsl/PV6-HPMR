function [ksp, fids, noise, rawKsp, sats] = ...
    shapeBrukerHpmr(ksp_in, method, idealHz, fMap, showImgs)
%SHAPEBRUKERHPMR shapes raw Bruker hpMR data into an MRI k-space
%
%   Usage: [ksp, fids, noise, rawKsp, sats] = ...
%               shapeBrukerHpmr(ksp_in, method, idealHz, fMap, showImgs)
%
%       where ksp_in is Bruker raw data reshaped to size
%               [nRead, nPhase, nRep] as by recoBrukerKspace
%             method is a Bruker method header structure
%               if omitted, readBrukerHeader('method') is called
%             idealHz contains the frequency offsets to use in IDEAL decomp
%               (only used if IDEAL echo shift is present)
%               if omitted, offsets are found in slice spectra, if present
%               if slice spectra aren't acquired, default PLHUAB are used
%             fMap is a 2D fieldmap matching image FOV, in units of Hz
%               (only used if IDEAL echo shift is present)
%             showImgs indicates plots to show:
%               if 1, display residual fids, slice spectra and data histograms
%               if 2, display residual fids and data histograms
%               default 1
%
%             ksp is reshaped k-space, IDEAL decomposed if appropriate
%             fids contains the pulse-acquire slice fids, if present
%             noise contains readouts acquired with negligible rf power
%               (raw, non-transformed noise, IDEAL decomposed if appropriate)
%             rawKsp is the raw k-space imaging readouts reshaped to size
%               [nRead, nSlice, nEcho, nRep]
%             sats contains the saturation echoes, if present
%
%   Literature:
%     Schulte RF, et al. MRM 2013
%     Wiesinger F, et al. MRM 2012
%
%   See also READBRUKERSTUDY, READBRUKERFID, RECOBRUKERKSPACE
%
%   07/2019, Keith Michel

funVersion = 'v20191106';

%% Check inputs
if nargin<1,	help(mfilename); return; end
if nargin<2,    readBrukerHeader('method'); end
if nargin<3,    idealHz = []; end
if nargin<4,    fMap = []; end
if nargin<5,    showImgs = 1; else, showImgs = round(double(showImgs(1))); end

if isdeployed
    fprintf(1, 'Entering %s (%s) at %s\n', mfilename, funVersion, datestr(now));
end

%% Parse scan type
if strcmpi(method.KAM_ExcMode, 'normal_exc')
    nSlice = method.PVM_SPackArrNSlices;
else
    if isfield(method, 'KAM_NSlices')
        nSlice = method.KAM_NSlices;
    else
        nSlice = 1;
    end
end

%% Read in appropriate design variables
if strcmpi(method.KAM_ReadMode, 'create2depi')
    % Readout designed on the fly using create2dEpi library
    try
        prefix   = evalin('caller', 'prefix;');
        tAcqPath = fullfile(prefix, 'tAcq');
    catch
        tAcqPath = 'tAcq';
    end
    fprintf(1, 'Attempting to generate create2dEpi recon parameters from tAcq file\n');
    if exist(tAcqPath, 'file')
        tAcqStruct = importdata(tAcqPath, '\n');
    else
        error('shapeBrukerHpmr:tAcqPath', 'Can''t find tAcq file at %s.', tAcqPath)
    end
    tAcq   = tAcqStruct.data(2:end).';
    idxAcq = tAcq > 0;
    tAcq   = tAcq(idxAcq);
    mtx    = method.KAM_EpiMatrix;
    pfy    = method.KAM_EpiPartialFourierY;
else
    % Readout pre-designed and saved in hpMR/gp
    hpmrPath = getenv('HPMR_PATH');
    if isempty(hpmrPath) && isdeployed
        if ispc
            [~, hpmrPath] = system('set PATH');
            hpmrPath = char(regexpi(hpmrPath, 'Path=(.*?);', 'tokens', 'once'));
        else
            hpmrPath = mfilename('fullpath');
        end
        while ~exist(fullfile(hpmrPath, 'hpMR', 'gp'), 'dir')
            tmp = fileparts(hpmrPath);
            if strcmp(tmp, hpmrPath)
                hpmrPath = [];
                break
            end
            hpmrPath = tmp;
        end
        if ~exist(fullfile(hpmrPath, 'hpMR', 'gp'), 'dir')
            hpmrPath = fullfile(getenv('HOME'), 'hpMR');
        end
    end
    fprintf(1, 'HPMR_PATH = %s\n', hpmrPath);

    readMode = [method.KAM_ReadMode, '.mat'];
    matPath  = fullfile(hpmrPath, 'hpMR', 'gp', readMode);
    if exist(matPath, 'file')
        load(fullfile(hpmrPath, 'hpMR', 'gp', readMode), ...
            'idxAcq', 'tAcq', 'mtx', 'pfy')
    else
        error('shapeBrukerHpmr:matPath', 'Can''t find mat file %s.', matPath)
    end
end
if ~exist('pfy', 'var'),    pfy = 1; end
mty  = mtx * pfy;
tAcq = tAcq - method.KAM_ReadEchoDelay; % Acq times relative to k-space center

%% Reshape data
if isfield(method, 'KAM_NEchoes')
    nEcho = method.KAM_NEchoes;
else
    nEcho = size(ksp_in, 2);
end
nRep = method.PVM_NRepetitions;

% Extract pulse acquire slice spectra if present
fids = [];
if strcmpi(method.KAM_FidEcho, 'yes')
    fids = squeeze(ksp_in(:,1,:));
    ksp_in(:,1,:) = [];
    if ~isfield(method, 'KAM_NEchoes')
        nEcho = nEcho - 1;
    end
end

% Extract saturation echo if present
sats  = [];
satOn = false;
if isfield(method, 'KAM_SatEcho')
    if strcmpi(method.KAM_SatEcho, 'yes')
        satOn = true;
        sats  = squeeze(ksp_in(:,1,:));
        ksp_in(:,1,:) = [];
        satFreq = [];
        if isfield(method, 'KAM_SatFreq')
            satFreq = method.KAM_SatFreq(1);
        end
    end
end
ksp = reshape(ksp_in, [], nSlice, nEcho, nRep);

% Extract noise echoes (with negligible excitation powers)
rfPow = method.KAM_VFAWatts;
noise = [];
if strcmpi(method.KAM_VFAOn, 'yes')
    if isfield(method, 'KAM_NAcqPoints')
        noise = ksp(1:method.KAM_NAcqPoints,:,rfPow<1e-3,:);
    else
        noise = ksp(:,:,rfPow<1e-3,:);
    end
    ksp(:,:,rfPow<1e-3,:) = [];
    nEcho = size(ksp, 3);
end

if method.KAM_IdealEchoShift > 1e-3
    idealOn = true;
    if isempty(idealHz)
        warning('shapeBrukerHpmr:idealNoOffsets', ...
            'hpMR scan has an IDEAL TE shift but no frequency offsets are provided. Decomposing PLHUAB.');
    end
else
    idealOn = false;
end
rawKsp = ksp;

%% Correct freq offsets using residual HP fid, if present
if isfield(method, 'KAM_NAcqPoints')
    nPoint = method.KAM_NAcqPoints;
else
    nPoint = method.PVM_EncMatrix(1);
end
if ~idealOn && nPoint > numel(idxAcq)+138
    % discard some points immediately after rewinders
    nJunk  = 10;
    resFid = ksp(numel(idxAcq)+nJunk:end,:,:,:);
    nFid   = size(resFid, 1);
    ksp(numel(idxAcq)+1:end,:,:,:) = [];
    ksp     = ksp(idxAcq,:,:,:);
    freqHz  = linspace(method.PVM_EffSWh/2, -method.PVM_EffSWh/2, 2048);
    if strcmpi(method.KAM_SpecExc, 'yes')
        txFreqs = method.KAM_VOFTxHertz;
    else
        txFreqs = zeros(1, nEcho);
    end
    txDegs  = method.KAM_VFADegrees;
    if numel(txFreqs) == 1
        txFreqs = txFreqs * ones(1, nSlice*nEcho);
    end
    if numel(txDegs) == 1
        txDegs = txDegs * ones(1, nSlice*nEcho);
    end
    if showImgs > 0
        figResFid = figure('Color', 'w', 'Name', 'hpMR Residual Spectra');
        axHandles = subplots(nSlice);
        legText   = {};
        cmap      = linspecer(nEcho+1, 'qualitative');
        if nEcho + 1 > 9
            whitebg(figResFid, 0.1*[1 1 1]);
        end
    end
    for ii = 1:nSlice
        for jj = 1:nEcho
            rewFid = zeros(2048, nRep);
            for kk = 1:nRep
                % Use reps with peak SNR >= 5
                curSpec = abs(ifftdim([resFid(:,ii,jj,kk); zeros(2048-nFid,1)], 1));
                if max(curSpec) >= 5*std(curSpec(1:128))/sqrt(2-pi/2)
                    rewFid(1:nFid,kk) = resFid(:,ii,jj,kk);
                end
            end
            % Sum over repetitions, 10 Hz line broadening
            rewFid  = sum(rewFid, 2);
            rewFid  = rewFid.' .* exp(-(0:2047)*10/method.PVM_EffSWh);
            rewSpec = abs(ifftdim(rewFid, 2));
            [~,idx] = max(rewSpec);
            offHz   = freqHz(idx);
            if abs(offHz) > 150
                fprintf(1, 'Slice %d, Echo %d (%+g Hz, %.02f Deg Tx), Residual relative spectrum offset %+g Hz (not applied)\n', ...
                    ii, jj, txFreqs(jj), txDegs(jj), offHz);
            else
                fprintf(1, 'Slice %d, Echo %d (%+g Hz, %.02f Deg Tx), Residual relative spectrum offset %+g Hz\n', ...
                    ii, jj, txFreqs(jj), txDegs(jj), offHz);
                ksp(:,ii,jj,:) = ksp(:,ii,jj,:) .* ...
                    repmat(exp(-2i*pi*tAcq'*offHz*1e-3), [1,1,1,nRep]);
            end
            if showImgs > 0
                figure(figResFid);
                semilogy(axHandles(ii) ,freqHz, rewSpec, 'Color', cmap(jj,:));
                legText{ii,jj} = sprintf('Echo %d: %+g Hz Offset, %.02f Deg Tx', jj, txFreqs(jj), txDegs(jj));
                hold on
            end
        end
    end
    ksp = reshape(ksp, mtx, mty, nSlice, nEcho, nRep);
    if satOn
        % discard 10 points immediately after rewinders, sum over repetitions
        resFid = sum(sats(numel(idxAcq)+nJunk:end,:), 2);
        sats(numel(idxAcq)+1:end,:) = [];
        sats     = sats(idxAcq,:);
        rewFid   = [resFid; zeros(2048-numel(resFid), 1)];
        rewFid   = rewFid.' .* exp(-(0:2047)*10/method.PVM_EffSWh);
        rewSpec  = abs(ifftdim(rewFid, 2));
        [mx,idx] = max(rewSpec);
        offHz    = freqHz(idx);
        rewSnr   = mx*sqrt(2-pi/2)/std(rewSpec(1:128));
        if rewSnr < 5 || abs(offHz) > 150
            fprintf(1, 'Saturation Echo (%+g Hz Tx), Residual relative spectrum offset %+g Hz (not applied)\n', ...
                satFreq, offHz);
        else
            fprintf(1, 'Saturation Echo (%+g Hz Tx), Residual relative spectrum offset %+g Hz\n', ...
                satFreq, offHz);
            sats = sats .* repmat(exp(-2i*pi*tAcq'*offHz*1e-3), [1,nRep]);
        end
        if showImgs > 0
            figure(figResFid);
            semilogy(axHandles, freqHz, rewSpec, 'Color', cmap(end,:));
            xlim(2e3*[-1 1]);
            xlabel('Freq (Hz)');
            set(gca, 'xdir', 'reverse')
            legText{ii,nEcho+1} = sprintf('Sat Echo, %+g Hz Tx Offset', satFreq);
            lHandle = legend(legText{ii,:});
            if ~verLessThan('matlab', '8.4')
                set(lHandle.BoxFace, 'ColorType', 'truecoloralpha', ...
                    'ColorData', uint8(255*[1;1;1;0.4]))
                grid on
            else
                set(lHandle, 'Color','none');
            end
            title(sprintf('Residual Spectra, Slice %d', ii), 'fontsize', 14)
        end
    elseif showImgs > 0
        figure(figResFid);
        xlim(2e3*[-1 1]);
        xlabel('Freq (Hz)');
        set(gca, 'xdir', 'reverse')
        lHandle = legend(legText);
        if ~verLessThan('matlab', '8.4')
            set(lHandle.BoxFace, 'ColorType', 'truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;0.4]))
            grid on
        else
            set(lHandle, 'Color','none');
        end
        title(sprintf('Residual Spectra, Slice %d', ii), 'fontsize', 14)
    end
elseif ~idealOn
    ksp(numel(idxAcq)+1:end,:,:,:) = [];
    ksp = reshape(ksp(idxAcq,:,:,:), mtx, mty, nSlice, nEcho, nRep);
else
    ksp(numel(idxAcq)+1:end,:,:,:) = [];
    ksp = ksp(idxAcq,:,:,:);
end


%% Process pulse-acquire spectra, if present
gB0 = method.PVM_FrqRef(1);  % MHz or Hz/ppm
if ~isempty(fids)
    f = linspace(method.PVM_EffSWh/2, -method.PVM_EffSWh/2, size(fids, 1));  % Hz
    t = (0:size(fids, 2)-1) * method.KAM_IdealRepDelay / 1e3; % s
    hpSpecs = ifftdim(fids, 1);
    subFreq = 4e3; % Hz
    idx     = ((f >= -subFreq/2) + (f <= subFreq/2)) > 1;
    ff      = f(idx);
    hpSpecs = abs(hpSpecs(idx,:));
    if showImgs == 1
        figure('Name', 'hpMR Slice Spectra')
        waterfall(ff, t, hpSpecs.');
        xlabel('Frequency (Hz)')
        ylabel('Time (sec)')
        set(gca, 'xdir', 'reverse', 'fontsize', 15)
        colormap(gca, brighten(flipud(cubehelix(64)), -0.8))
    end
    if idealOn && isempty(idealHz)
        % IDEAL recon using shifts from hpFid spectra
        hpSpecs = flipud(sum(hpSpecs, 2));
        ff      = fliplr(ff);
        ppm     = method.PVM_FrqWorkPpm(1) + ff/gB0;

        % The following works in LDH phantoms,
        % but may not be suitable for in vivo studies
        [~, pk] = findpeaks(hpSpecs, ...
            'NPeaks',	4, ...      % pyr, lac, hyd, urea
            'SortStr',	'descend');
        pk       = ppm(pk);
        ppmExp   = [171, 183.2, 179, 163.5];
        idealPpm = zeros(1, numel(ppmExp));
        for ii = 1:numel(ppmExp)
            [~,idx]      = min(abs(pk-ppmExp(ii)));
            if abs(pk(idx)-ppmExp(ii)) < 1.5
                idealPpm(ii) = pk(idx);
            else
                idealPpm(ii) = ppmExp(ii);
            end
            pk(idx) = [];
        end
        idealPpm = [idealPpm, 177, 161]; % ala, bic
    end
else
    if idealOn && isempty(idealHz)
        warning('shapeBrukerHpmr:idealNoOffsetsNoFids', ...
            'hpMR IDEAL scan acquired without slice spectra, using default ppm offsets for PLHUAB decomp.')
        idealPpm = [171, 183.2, 179.5, 163.5, 177, 161]; % PLHUAB
    end
end

%% Perform IDEAL decomposition
if idealOn
    if isempty(idealHz)
        idealHz = (idealPpm - method.PVM_FrqWorkPpm(1)) * gB0;
    elseif ~exist('idealPpm', 'var')
        idealPpm = idealHz / gB0 + method.PVM_FrqWorkPpm(1);
    end
    fprintf(1, 'Species: Pyr     Lac     Hyd     Urea    Ala     BiC\n');
    fprintf(1, '    PPM: %-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f\n', idealPpm);
    fprintf(1, '     Hz: %-8.1f%-8.1f%-8.1f%-8.1f%-8.1f%-8.1f\n\n', idealHz);
    tes = method.PVM_EchoTime + (0:nEcho-1)*method.KAM_IdealEchoShift;
    % Fieldmap demodulation
    if ~isempty(fMap)
        fprintf(1, 'Demodulating fieldmap from individual images\n');
        pMap = zeros(mtx, mty, nEcho);
        for ii = 1:nEcho
            pMap(:,:,ii) = exp(-1i * 2e-3 * pi * tes(ii) * fMap);
        end
        pMap = repmat(reshape(pMap, mtx, mty, 1, nEcho), 1, 1, nSlice, 1);
        tImg = ifftdim(reshape(ksp, mtx, mty, nSlice, nEcho, nRep), 1:2);
        tImg = tImg .* repmat(pMap, 1, 1, 1, 1, nRep);
        ksp  = reshape(fftdim(tImg, 1:2), mtx*mty, nSlice, nEcho, nRep);
    end
    % IDEAL decomp of kspace data w dwell time and t2* corrections
    t2s   = 20; % ms
    nChem = numel(idealHz);
    kd    = zeros(mtx, mty, nSlice, nChem, nRep);
    for ii = 1:nRep
        for jj = 1:nSlice
            curkd = idDecomp(squeeze(ksp(:,jj,:,ii)), ...
                idealHz, tes, tAcq, t2s, false);
            kd(:,:,jj,:,ii) = reshape(curkd, mtx, mty, 1, nChem, 1);
        end
    end
    ksp = kd;
    % Attempt IDEAL decomp of acquired noise readouts
    try
        noise = reshape(noise(idxAcq,:,:,:), mtx*mty, []);
        nNoi  = floor(size(noise, 2) / nEcho);
        noise = reshape(noise(:,1:nEcho*nNoi), mtx*mty, nEcho, nNoi);
    catch
        warning('shapeBrukerHpmr:idealNoise', ...
            'Failed to reshape noise for IDEAL decomposition, discarding noise')
        noise = [];
    end
    if ~isempty(noise)
        tmp = zeros(mtx, mty, nChem, nNoi);
        for ii = 1:nNoi
            curkd = idDecomp(squeeze(noise(:,:,ii)), ...
                idealHz, tes, tAcq, t2s, false);
            tmp(:,:,:,ii) = reshape(curkd, mtx, mty, [], 1);
        end
        noise = tmp;
    end
    nEcho = nChem;
end

%% Image phase corrections
% Apply k-space shift to phase encode dimension
kShift = -method.PVM_Phase1Offset ./ method.PVM_Fov(2);
for ii = 1:nSlice
    ksp(:,:,ii,:,:) = shiftKspace(ksp(:,:,ii,:,:), [0 kShift]);
end
% Zero-pad partial Fourier data along phase encode dimension
nPfy = mtx * (1-pfy);
ksp  = [zeros(mtx,nPfy,nSlice,nEcho,nRep), ksp];
% Apply N/2 ghost correction for symmetric EPI readouts
ghostCorr = false;
switch lower(method.KAM_ReadMode(1:6))
    case 'symepi'
        ghostCorr = true;
    case 'create'
        if strcmpi(method.KAM_EpiType, 'symmetric')
            ghostCorr = true;
        end
end
if ghostCorr
    ksp = ghostCorrEpi(ksp);
end

%% Show histograms of noise and k-space signal, if noise is acquired
if ~isempty(noise) && showImgs > 0
    nBins = 30;
    figure('Color', 'w', 'Name', 'hpMR Noise', 'Position', [50 50 900 600])
    subplot 121
    if verLessThan('matlab', '8.4')
        [n,x] = hist(real(noise(:)), nBins);
        b = bar(x, n, 'hist');
        set(b, 'FaceColor', [0 0.447 0.741], 'FaceAlpha', 0.5);
    else
        histogram(real(noise), nBins);
    end
    hold on
    if verLessThan('matlab', '8.4')
        [n,x] = hist(imag(noise(:)), nBins);
        b = bar(x, n, 'hist');
        set(b, 'FaceColor', [0.85 0.325 0.098], 'FaceAlpha', 0.5);
    else
        histogram(imag(noise), nBins);
    end
    xlabel('Real / Imaginary Value')
    title(sprintf('Noise (Inverse CV = %.04f)', ...
        mean(abs(noise(:)))/std(abs(noise(:)))))
    set(gca, 'FontSize', 15)
    if ~verLessThan('matlab', '8.4')
        set(gca, 'YScale', 'log')
    end
    subplot 122
    if verLessThan('matlab', '8.4')
        [n,x] = hist(real(ksp(:)), nBins);
        b = bar(x, n, 'hist');
        set(b, 'FaceColor', [0 0.447 0.741], 'FaceAlpha', 0.5);
    else
        histogram(real(ksp), nBins);
    end
    hold on
    if verLessThan('matlab', '8.4')
        [n,x] = hist(imag(ksp(:)), nBins);
        b = bar(x, n, 'hist');
        set(b, 'FaceColor', [0.85 0.325 0.098], 'FaceAlpha', 0.5);
    else
        histogram(imag(ksp), nBins);
    end
    title('K-space')
    xlabel('Real / Imaginary Value')
    set(gca, 'FontSize', 15)
    if ~verLessThan('matlab', '8.4')
        set(gca, 'YScale', 'log')
    end
    movegui(gcf, 'south')
end

if isdeployed
    fprintf(1, 'Leaving %s at %s\n', mfilename, datestr(now));
end

end
