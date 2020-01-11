function [img, ksp, hpFids, hpNoise, hpKraw, hpSats] = ...
        recoBrukerKspace(brukerFid, acqp, method, showImgs, idealHz)
%RECOBRUKERKSPACE shapes raw Bruker MR data into an MRI k-space and applies
%   inverse discrete Fourier transform where appropriate
%
%   Usage: 
%     Generic MR Experiment:
%       [img, ksp] = recoBrukerKspace(brukerFid, acqp, method, showImgs)
%     hpMR Experiment:
%       [img, ksp, hpFids, hpNoise, hpKraw, hpSats] = ...
%           recoBrukerKspace(brukerFid, acqp, method, showImgs, idealHz)
%
%       where brukerFid is Bruker raw data
%               if scalar or string, readBrukerFid is called in that subdir
%               if omitted, readBrukerFid is called in working directory
%             acqp is a Bruker acqp header structure
%               if omitted, readBrukerHeader('acqp') is called
%             method is a Bruker method header structure
%               if omitted, readBrukerHeader('method') is called
%             showImgs is a logical, if true display transformed data
%               default true 
%
%       for hpMR experiments:
%             showImgs indicates plots to show:
%               if 1, display all plots, ifdeployed save outputs to .mat file
%               if 2, display only residual fids and data histograms
%               default 1 (passed directly to SHAPEBRUKERHPMR)
%             idealHz contains the frequency offsets to use in IDEAL decomp
%               default [] (PLHUAB)
%
%   See also READBRUKERHEADER, READBRUKERFID, SHAPEBRUKERHPMR
%
%   07/2019, Keith Michel

funVersion = 'v20191103';

%% Check inputs
if ~nargin
    if isdeployed, brukerFid = uigetdir('', 'Select Bruker Experiment Directory');
    else,          brukerFid = readBrukerFid; end
end
if isscalar(brukerFid) || ischar(brukerFid)
    prefix    = num2str(brukerFid);
    brukerFid = readBrukerFid(brukerFid);
else
    prefix = [];
end
if nargin<2,         acqp     = []; end
if isempty(acqp),    acqp     = readBrukerHeader(fullfile(prefix, 'acqp')); end
if nargin<3,         method   = []; end
if isempty(method),  method   = readBrukerHeader(fullfile(prefix,'method')); end
if nargin<4,         showImgs = 1; end
if nargin<5,         idealHz  = []; end

if isdeployed
    if showImgs == 1
        diary(fullfile(prefix, sprintf('%s-log_%s.txt', mfilename, datestr(now, 30))));
    end
    fprintf(1, 'Entering %s (%s) at %s\n', mfilename, funVersion, datestr(now));
end

hpFids  = [];
hpNoise = [];
hpKraw  = [];
hpSats  = [];

%% Parse data size
dim    = acqp.ACQ_dim;
nPoint = prod(acqp.ACQ_size)/2;
nEcho  = acqp.NECHOES;
nImg   = acqp.NI;
nRep   = acqp.NR;
nSlice = acqp.NSLICES;
if any(strcmpi(acqp.ACQ_dim_desc, 'spatial'))
    nPhase  = acqp.ACQ_size(2); 
    nPoint = nPoint / nPhase;
    if dim == 3
        nSlice = acqp.ACQ_size(3);
        nPoint = nPoint / nSlice;
        nImg   = nImg * nSlice;
    end
else
    nPhase = 0;
end

nRead = 128 * ceil(nPoint/128);

%% Custom reordering routine for individual pulseprograms
ppg = lower(acqp.PULPROG);
if strcmp(ppg, 'hpmr.ppg')
    if strcmp(method.KAM_ReadMode, 'Normal_Read')
        ppg = 'flash.ppg';
    end
elseif strcmp(ppg, 'spflash.ppg')
    nRep = nPhase * nRep;
    dim  = 1;
end
% TODO: csi
switch ppg
    % General imaging methods
    case {'flash.ppg', 'fspgr.ppg', 'vfaflash.ppg', 'fspgrvfa.ppg', ...
            'mge.ppg', 'rare.ppg', 'rarevtr.ppg', 'msme.ppg', ...
            'fcflash.ppg', 'fcflashangio.ppg'}
        % Echoes must be for separate phase encodes or separate images
        if nImg == nSlice
            nPhase = nPhase/nEcho;
        elseif nImg ~= nSlice * nEcho
                error('recoBrukerKspace:nonUniformEchoAveraging', ...
                    'Multiecho scans must iterate phase encodes (RARE) or echo images (MSME, MGE).')
        end
        % Try initial reshape
        try
            ksp = reshape(brukerFid, ...
                nRead, nEcho, nSlice, nPhase, nRep);
        catch
            error('recoBrukerKspace:flashRareEtc', ...
                'Failed to reshape %s MR data', ppg)
        end
        [~,ordObj] = sort(acqp.ACQ_obj_order);
        if dim == 2 && nImg == nSlice
            % One echo image, ACQ_obj_order is for slices
            ksp = ksp(:,:,ordObj,:,:);
        end
        if nImg == nSlice * nEcho && nEcho > 1
            % Separate echo images, ACQ_obj_order is for echoes and slices
            if dim == 2
                ksp = reshape(ksp, ...
                    nRead, nEcho*nSlice, nPhase, nRep);
                ksp = ksp(:,ordObj,:,:);
                ksp = reshape(ksp, ...
                    nRead, nEcho, nSlice, nPhase, nRep);
            else
                ksp = ksp(:,ordObj,:,:,:);
            end
            ksp = permute(ksp, [1, 4, 3, 2, 5]);
        else
            % Phase encoded echoes
            if dim == 2
                ksp = reshape(permute(ksp, [1, 2, 4, 3, 5]), ...
                    nRead, nEcho*nPhase, nSlice, nRep);
            else
                ksp = reshape(ksp, nRead, nEcho*nPhase, nSlice, nRep);
            end
        end
        [~,ordPhase] = sort(acqp.ACQ_spatial_phase_1);
        ksp = ksp(:,ordPhase,:,:);
        if dim == 3
            [~,ordPhase] = sort(acqp.ACQ_spatial_phase_2);
            ksp = ksp(:,:,ordPhase,:);
        end
        % Apply k-space shift to phase-encode dimension
        kShift = -acqp.ACQ_phase1_offset ./ acqp.ACQ_fov(2) / 10;
        for ii = 1:nSlice
            ksp(:,:,ii,:,:) = shiftKspace(ksp(:,:,ii,:,:), [0, kShift(ii)]);
        end
        
    % Simple spectroscopy methods
    case {'spflash.ppg', 'singlepulse.ppg', 'powcalsinglepulse.ppg', ...
            'nspect.ppg', 'press.ppg'}
        % Try initial reshape
        try
            ksp = reshape(brukerFid, nRead, nRep);
        catch
            error('recoBrukerKspace:nspectPressEtc', ...
                'Failed to reshape %s MR data', ppg)
        end
        
    % Bruker hpMR
    case 'hpmr.ppg'
        % Try initial reshape
        try
            ksp = reshape(brukerFid, nRead, nPhase, nRep);
        catch
            error('recoBrukerKspace:hpMR', ...
                'Failed to reshape %s MR data', ppg)
        end
        [ksp, hpFids, hpNoise, hpKraw, hpSats] = shapeBrukerHpmr(...
            ksp, method, idealHz, [], showImgs);
        
    % Unsupported methods
    otherwise
        error('recoBrukerKspace:ppgNotSupported', ...
            'This function currently doesn''t support pulseprogram %s', ppg)
end

%% Remove zero-padding along readout
if mod(nPoint, 128)
    ksp(acqp.ACQ_size(1)/2+1:end,:,:,:) = [];
end

%% Reconstruct and show images
if ~isempty(hpFids) && showImgs == 1
    plotSpec(ifftdim(hpFids, 1), acqp.SW_h, true);
    set(gcf, 'Name', sprintf('hpMR Slice Spectra - %s', acqp.ACQ_scan_name))
end
img = ifftdim(ksp, 1:dim);
if showImgs == 1
    if dim == 1
        plotSpec(img, acqp.SW_h);
    else
        showMontage(img);
    end
    set(gcf, 'Name', sprintf('%s Images - %s', mfilename, acqp.ACQ_scan_name))
    movegui(gcf, 'center')
end

if isdeployed
    if showImgs == 1
        % save(fullfile(prefix, sprintf('%s-workspace_%s', ...
        %     mfilename, datestr(now, 30))), 'img', 'ksp', 'hp*')
        if strcmp(ppg, 'hpmr.ppg')
            saveVars = {'img', 'ksp', 'funVersion', 'acqp', 'method'...
            'hpFids', 'hpNoise', 'hpKraw', 'hpSats', 'idealHz'};
        else
            saveVars = {'img', 'ksp', 'funVersion', 'acqp', 'method'};
        end
        uisave(saveVars, fullfile(prefix, sprintf('%s-workspace_%s_%s', ...
            mfilename, strtok(ppg, '.'), datestr(now, 30))));
        if dim > 1
            viewover(squeeze(img), 'title', acqp.ACQ_scan_name)
        end
        waitfor(gcf)
    end
    fprintf(1, 'Leaving %s at %s\n', mfilename, datestr(now));
end


end
