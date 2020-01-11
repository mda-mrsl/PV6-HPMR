function outStruct = hpFidEcho(inDir,params)
%HPFIDECHO reads in and analyze slice spectra from hpMR experiment
%
%   Usage: outStruct = hpFidEcho(inDir, params)
%
%   Inputs:
%     inDir  - Integer or char of IDEAL scan directory
%     params - Structure outlining analysis parameters with elements:
%       searchBands - Array of spectral bands to find peaks (ppm)
%         default = [169, 173; 182, 185] (pyr; lac)
%       lbHz        - Exponential line broadening factor (Hz)
%         default = 10 
%       qMethod     - String specifying quantification method
%         default = 'fwhm'    - FWHM integration
%                   'genspec' - Fit peak amplitudes in freq domain
%       noiseRegion - 2-element vector specifying limits of noise
%                     region used for baseline correction (ppm)
%         default = [195 200]
%       timeRegion  - 2-element vector specifying time region for analysis (sec)
%         default = [0 150]
%       ppmRange    - 2-element vector specifying spectral coordinates for 
%                     cropping waterfall plot
%         default = [160 190]
%       showPlots   - Logical, Show result plots
%         default = true
%
%   Outputs: Structure containing elements:
%     spectra     - Complex slice spectra over time
%     avgSpectrum - Time-averaged magnitude spectrum
%     ppmAxis     - Spectral axis coordinates (ppm)
%     timeAxis    - Time axis coordinates (sec)
%     timeCurves  - Array of individual species time curves (normalized)
%     snrCurves   - Array of individual species SNR curves
%     ratioCurves - Array of cumulative ratiometric time curves for each
%                   species, first species is assumed to be precursor
%     chemShifts  - Vector of individual species coordinates (ppm)
%     chemShiftsPw- Array of point-wise individual species coordinates
%                   over time (ppm)
%     lineWidths  - Vector of individual species linewidths (ppm)
%     lineWidthsPw- Array of point-wise individual species linewidths
%                   over time (ppm)
%     spSnr       - SNR of spectrum (peak SNR for each slice)
%     headers     - Combined acqp and method headers
%     funVersion  - faFidEcho function version
%
%   Literature: 
%
% 02/2018, Keith Michel

%if ~nargin, help(mfilename); return; end;

funVersion = 'v20190910';
if isdeployed
    fprintf(1, 'Entering %s (%s) at %s\n', mfilename, funVersion, datestr(now));
end

%% Parse inputs
if ~nargin, inDir = uigetdir(); end
if ~exist('params','var'),      params = []; end
parsein(params);

if isscalar(inDir),             inDir = num2str(inDir); end
if ~exist(inDir,'dir')
    error('Specified directory does not exist');
end

nShifts = size(searchBands,1);
if size(searchBands,2) ~= 2
    error('Chemical shift search bands must be an Nx2 array');
end

if numel(noiseRegion) > 2
    warning('faFidEcho:noiseRegion','numel(noiseRegion) = %d > 2. Using first two values.',numel(noiseRegion));
    noiseRegion = noiseRegion(1:2)
elseif numel(noiseRegion) < 2
    error('Noise region must contain 2 elements.')
end

if numel(timeRegion) > 2
    warning('faFidEcho:timeRegion','numel(timeRegion) = %d > 2. Using first two values.',numel(timeRegion));
    timeRegion = timeRegion(1:2)
elseif numel(timeRegion) < 2
    error('Time region must contain 2 elements.')
end

if numel(ppmRange) > 2
    warning('faFidEcho:ppmRange','numel(ppmRange) = %d > 2. Using first two values.',numel(ppmRange));
    ppmRange = ppmRange(1:2)
elseif numel(ppmRange) < 2
    error('Ppm range must contain 2 elements.')
end


%% Read headers
headers = struct();
method = readBrukerHeader(fullfile(inDir,'method'));
acqp = readBrukerHeader(fullfile(inDir,'acqp'));
if isempty(strfind(acqp.ACQ_method,'hpMR'))
    error('Not a hpMR experiment')
end
% Fill with method info
tmpNames = fieldnames(method);
for ii = 1:numel(tmpNames)
    headers.(tmpNames{ii}) = method.(tmpNames{ii});
end
if strcmp(method.KAM_FidEcho,'No')
    error('Fid echo not acquired in specified scan');
end
% Fill with acqp info
tmpNames = fieldnames(acqp);
duplicateNames = {'TITLE','JCAMPDX','DATATYPE','ORIGIN','OWNER'};
for ii = 1:numel(duplicateNames)
    if(isfield(acqp,duplicateNames{ii}))
        tmpNames(strcmp(tmpNames,duplicateNames{ii})) = [];
    end
end
for ii = 1:numel(tmpNames)
    if(isfield(headers,tmpNames{ii}))
        warning('*WARNING* %s is a field name in both the method and acqp file. Using the value from the acqp file/n',tmpNames{ii})
    end
    headers.(tmpNames{ii}) = acqp.(tmpNames{ii});
end

clear acqp method

%% Read and reshape data
nReps = headers.PVM_NRepetitions;
nSlices = headers.NSLICES;
timeAxis = 0:headers.KAM_IdealRepDelay:headers.KAM_IdealRepDelay*(nReps-1);
timeAxis = timeAxis/1000;   % sec
if timeRegion(end) > timeAxis(end)
    warning('faFidEcho:timeRegion','timeRegion(end) = %g, but scan is only %g seconds long',timeRegion(2),timeAxis(end));
    timeRegion(end) = timeAxis(end)
end
timeRange(1) = find( timeAxis >= timeRegion(1), 1);
timeRange(2) = find( timeAxis >= timeRegion(2), 1);
nPtsRaw = headers.ACQ_size(1)/2;
% specRes = 10000/512;    % desired spectral resolution, Hz/pt
% specRes = 31250/625;
specRes = 1e6/48/2048;
nPts = round( headers.SW_h / specRes );
f = fopen(fullfile(inDir,'fid'));
dat = fread(f,inf,'int32');
fclose(f);
fid = dat(1:2:end) + 1i*dat(2:2:end);
fidr = reshape(fid,[],nSlices,headers.ACQ_size(2),nReps);
spFids = squeeze( permute( fidr(1:nPtsRaw,:,1,:), [1 3 4 2]));
clear dat fid fidr

%% Line broadening and IFT
spFids = spFids .* repmat( exp( -lbHz * (0:1:nPtsRaw-1)/headers.SW_h).', [1, nReps, nSlices]);
if nPts > nPtsRaw, spFids(end+1:nPts,:,:) = 0; end
spectra = flipdim(ifftdim(spFids(1:nPts,:,:), 1), 1);

%% Handle spectral axis
gB0 = headers.PVM_FrqRef(1);  % MHz or Hz/ppm
ppmCenter = 1e6*(headers.PVM_FrqWork(1) - headers.PVM_FrqRef(1))/gB0;
ppmAxis = ppmCenter + linspace(-headers.SW/2,headers.SW/2,nPts);
if any(searchBands(:)<ppmAxis(1)) || any(searchBands(:)>ppmAxis(end))
    error('Chemical shift search bands must be within scan bandwidth [%.1f to %.1f]',ppmAxis(1),ppmAxis(end));
end
idxNoise = [find( ppmAxis >= noiseRegion(1), 1), ...
    find( ppmAxis >= noiseRegion(2), 1)];
idxSearch = zeros(nShifts,2);
for ss = 1:nShifts
    idxSearch(ss,1) = find( ppmAxis >= searchBands(ss,1), 1);
    idxSearch(ss,2) = find( ppmAxis >= searchBands(ss,2), 1);
end

%% Iterate through repetitions, get individual species time curves
timeCurves   = zeros(nShifts,nReps,nSlices);
snrCurves    = zeros(nShifts,nReps,nSlices);
chemShiftsPw = zeros(nShifts,nReps,nSlices);
lineWidthsPw = zeros(nShifts,nReps,nSlices);
ratioCurves  = zeros(nShifts-1,nReps,nSlices);   % first shift is precursor (e.g. pyr)
corrSpectra  = zeros(size(spectra));
noise        = abs(spectra(idxNoise(1):idxNoise(2),:,:));
stdNoise     = std(noise(:))/sqrt(2-pi/2);
switch qMethod
    case 'fwhm' % FWHM integration of magnitude peak in each search band
        for sl = 1:nSlices
            for tt = 1:nReps
                curSpectrum = abs(spectra(:,tt,sl));
                for ss = 1:nShifts
                    [timeCurves(ss,tt,sl), chemShiftsPw(ss,tt,sl), lineWidthsPw(ss,tt,sl)] = ...
                        fwhm( curSpectrum, idxSearch(ss,:) );
                    snrCurves(ss,tt,sl)    = timeCurves(ss,tt,sl) / ...
                        lineWidthsPw(ss,tt,sl) / stdNoise;
                    chemShiftsPw(ss,tt,sl) = ppmAxis( chemShiftsPw(ss,tt,sl));
                    lineWidthsPw(ss,tt,sl) = lineWidthsPw(ss,tt,sl) * specRes / gB0;
                end
                for ss = 2:nShifts
                    ratioCurves(ss-1,tt,sl) = sum( timeCurves(ss,1:tt,sl) ) / ...
                        ( sum( timeCurves(ss,1:tt,sl) ) + sum( timeCurves(1,1:tt,sl) ) );
                end
                corrSpectra(:,tt,sl) = curSpectrum;
            end
        end
    case 'genspec'
        opts = optimset('lsqcurvefit');
        opts = optimset(opts,'MaxIter',5e4,'MaxFunEvals',5e4,...
            'TolFun',1e-5,'TolX',1e-5,'display','off');
%         tic
        for sl = 1:nSlices
            for tt = 1:nReps
                curNoise = mean( abs( spectra(idxNoise(1):idxNoise(2),tt,sl) ) );
                curSpectrum = abs(spectra(:,tt,sl)) - curNoise;
                guessPars = repmat([0 0 20/specRes].',1,nShifts);
                ubPars = repmat([Inf 0 100/specRes].',1,nShifts);
                lbPars = repmat([0 0 .2/specRes].',1,nShifts);
                for ss = 1:nShifts
                    [pk,loc] = max(curSpectrum(idxSearch(ss,1):idxSearch(ss,2)));
                    guessPars(1,ss) = pk;
                    guessPars(2,ss) = loc + idxSearch(ss,1) - 1;
                    ubPars(2,ss) = idxSearch(ss,2);
                    lbPars(2,ss) = idxSearch(ss,1);
                end
                fitPars = lsqcurvefit( @genspec, guessPars, ppmAxis, ...
                    curSpectrum, lbPars, ubPars, opts);
                fitPars = reshape( fitPars, 3, nShifts);
%                 figure(9); clf; plot(ppmAxis,curSpectrum); hold on; plot(ppmAxis,genspec(fitPars,ppmAxis)); pause(0.2);
                corrSpectra(:,tt) = genspec(fitPars,ppmAxis);
                for ss = 1:nShifts
                    [timeCurves(ss,tt,sl), chemShiftsPw(ss,tt,sl), lineWidthsPw(ss,tt,sl)] = ...
                        fwhm( genspec(fitPars(:,ss),ppmAxis), idxSearch(ss,:) );
                    chemShiftsPw(ss,tt,sl) = ppmAxis( chemShiftsPw(ss,tt,sl));
                    lineWidthsPw(ss,tt,sl) = lineWidthsPw(ss,tt,sl) * specRes / gB0;
                end
                for ss = 2:nShifts
                    ratioCurves(ss-1,tt,sl) = sum( timeCurves(ss,1:tt,sl) ) / ...
                        ( sum( timeCurves(ss,1:tt,sl) ) + sum( timeCurves(1,1:tt,sl) ) );
                end
            end
        end
%         toc
%         for tt = 1:nReps
%             for ss = 2:nShifts
%                 ratioCurves(ss-1,tt) = sum( timeCurves(ss,1:tt) ) / ...
%                     ( sum( timeCurves(ss,1:tt) ) + sum( timeCurves(1,1:tt) ) );
%             end
%         end
end

%% Get average spectrum, chem shifts and normalized species time curves
subTimeRange = timeRange(1):timeRange(2);
avgSpectrum = zeros(size(spectra,1),nSlices);
chemShifts = zeros(nShifts,nSlices);
lineWidths = zeros(nShifts,nSlices);
spSnr = zeros(1,nSlices);
fString = [];
for sl = 1:nSlices
    avgSpectrum(:,sl) = sum( abs(corrSpectra(:,subTimeRange,sl)), 2);
    avgSpectrum(:,sl) = avgSpectrum(:,sl) - mean(avgSpectrum(idxNoise(1):idxNoise(2),sl));
    % avgSig = 0;
    for ss = 1:nShifts
        % pk = max(avgSpectrum(idxSearch(ss,1):idxSearch(ss,2)));
        % avgSig = avgSig + pk;
        [~,chemShifts(ss,sl),lineWidths(ss,sl)] = fwhm( avgSpectrum, idxSearch(ss,:));
        chemShifts(ss,sl) = ppmAxis( chemShifts(ss,sl));
        lineWidths(ss,sl) = lineWidths(ss,sl) * specRes / gB0;
    end
    % spSnr(sl) = avgSig / stdNoise;
    maxSig = max( max( abs(corrSpectra(idxSearch(1,1):idxSearch(1,2),:,sl))));
    spSnr(sl) = maxSig / stdNoise;      % peak SNR for precursor over time
    timeCurves(:,:,sl) = timeCurves(:,:,sl) / ...
        max( max(timeCurves(:,:,sl)));  % normalize time curves to max signal
    fString = [fString ...
        sprintf('Slice %d: nLac %g, peak pyr SNR %g, pyr linewidth %g ppm (%g ms T2*)\n',...
        sl, ratioCurves(1,subTimeRange(end),sl),spSnr(sl), ...
        lineWidths(1,sl), 1e3/gB0/pi/lineWidths(1,sl))];
end
timeCurves   = timeCurves(:,subTimeRange,:);
ratioCurves  = ratioCurves(:,subTimeRange,:);
chemShiftsPw = chemShiftsPw(:,subTimeRange,:);
lineWidthsPw = lineWidthsPw(:,subTimeRange,:);

%% Plot results
if showPlots
    
    if nShifts == 2
        colors = [0 0.447 0.741; 0.85 0.325 0.098];
    elseif ~verLessThan('matlab', '8.5')
        colors = lines(nShifts);
    elseif exist('brewermap','file')
        colors = brewermap(nShifts, 'Dark2');

    end
    
    figResults = struct;
    figWaterfall = struct;
    
    subPpmRange = find(ppmAxis>=ppmRange(1),1) : find(ppmAxis>=ppmRange(2),1);
    subPpmAxis = ppmAxis(subPpmRange);

    for sl = 1:nSlices
       
        curFigResults = feval('figure', 'Name', sprintf( ...
            'Slice Spectra Results - Slice %d - nLac %0.2f - Peak SNR %0.2f', ...
            sl, ratioCurves(1,end,sl), spSnr(sl)), 'Position',...
            [50 50 1800 600], 'Color', 'w');

        % time-spectra image
        subplot 131
        curSpectra = abs(spectra(:,:,sl));
        imagesc( timeAxis, ppmAxis, squeeze(curSpectra) );
        caxis([0 max( curSpectra(:))/3]);
        if exist('cubehelix','file')
            cmap = cubehelix(64);
        else
            cmap = bone(64);
        end
        colormap(gca,cmap);
        ylim([subPpmAxis(1) subPpmAxis(end)])
        set(gca,'ydir','reverse','fontsize',14);
        ylabel('Chemical Shift (ppm)');
        xlabel('Time (sec)');
        title('Spectra','fontsize',16);
        hold on
        line(timeAxis(timeRange(1)*[1 1]),[ppmAxis(1) ppmAxis(end)],...
            'color',[1 1 1],'linewidth',1.6);
        line(timeAxis(timeRange(2)*[1 1]),[ppmAxis(1) ppmAxis(end)],...
            'color',[1 1 1],'linewidth',1.6);
        for ss = 1:nShifts
            line([timeAxis(1) timeAxis(end)],ppmAxis(idxSearch(ss,1))*[1 1],...
                'color',colors(ss,:),'linewidth',1.6);
            line([timeAxis(1) timeAxis(end)],ppmAxis(idxSearch(ss,2))*[1 1],...
                'color',colors(ss,:),'linewidth',1.6);
        end

        % time-averaged spectrum
        subplot 132
        plot(ppmAxis,avgSpectrum(:,sl),'color',0.1*ones(1,3));
        hold on
        for ss = 1:nShifts
            curBand = idxSearch(ss,1):idxSearch(ss,2);
            area(ppmAxis(curBand),avgSpectrum(curBand,sl),...
                'facecolor',colors(ss,:),'edgecolor',colors(ss,:));
        end
        area(ppmAxis(idxNoise(1):idxNoise(2)),avgSpectrum(idxNoise(1):idxNoise(2),sl),...
            'facecolor',[.7 .7 .7],'edgecolor',[.7 .7 .7]);
        xlim([subPpmAxis(1) subPpmAxis(end)])
        set(gca,'xdir','reverse','fontsize',14); grid on
        xlabel('Chemical Shift (ppm)');
        ylabel('Mean Signal (AU)');
        title('Time-averaged Spectrum','fontsize',16);

        % species time curves
        subplot 133
        for ss = 1:nShifts
            plot(timeAxis(subTimeRange),timeCurves(ss,:,sl),...
                'color',colors(ss,:),'linewidth',2.2);
            hold on
        end
        for ss = 2:nShifts
            plot(timeAxis(subTimeRange),ratioCurves(ss-1,:,sl),...
                'color',colors(ss,:),'linewidth',2.2,'linestyle','--');
        end
        if numel(subTimeRange) > 1
            xlim(timeAxis([subTimeRange(1) subTimeRange(end)]))
        end
        ylim([-0.1 1.1])
        set(gca,'fontsize',14); grid on
        xlabel('Time (sec)')
        ylabel('Relative Signal')
        title('Species Time Curves','fontsize',16);
        
        if exist('export_fig', 'file')
            export_fig(curFigResults, fullfile(inDir, sprintf( ...
                'faFidEcho_resultsPlot-Slice%d',sl)),'-a1');
        else
            print(curFigResults, fullfile(inDir, sprintf( ...
                'faFidEcho_resultsPlot-Slice%d',sl)),'-dpng');
        end
        eval(sprintf('figResults.slice%d = curFigResults;',sl))

        % waterfall plot
        if showPlots == 1
            
            curFigWaterfall = feval('figure', 'Name',...
                sprintf('Waterfall Plot - Slice %d',sl), 'Position',...
                [50 50 800 800], 'Color', 'w');
            
            wfsubTimeRange = subTimeRange(1:2:end); % less cluttered plot
            waterfalls = corrSpectra(subPpmRange,wfsubTimeRange,sl);
            waterfalls = waterfalls - min(waterfalls(:));
            waterfalls = waterfalls / max(waterfalls(:));
            normAvgSpectrum = avgSpectrum(subPpmRange,sl);
            normAvgSpectrum = normAvgSpectrum - min(normAvgSpectrum);
            normAvgSpectrum = normAvgSpectrum / max(normAvgSpectrum);
            waterfall(subPpmAxis,timeAxis(wfsubTimeRange),waterfalls.');
            set(gca,'xdir','reverse','fontsize',14); grid on
            if exist('cubehelix','file')
                cmap = brighten( flipud( cubehelix(64)),-.8);
            elseif exist('brewermap','file')
                cmap = brighten( brewermap(64,'Blues'),-.85);
            else
                cmap = flipud( winter(64));
            end
            colormap(gca,cmap);
            hold on
            plot3(subPpmAxis,timeAxis(wfsubTimeRange(end))*ones(1,numel(subPpmAxis)),...
                normAvgSpectrum,'color',[.2 .2 .2],'linewidth',2);
            for ss = 1:nShifts
                plot3(ppmRange(1)*ones(1,numel(wfsubTimeRange)),timeAxis(wfsubTimeRange),...
                    timeCurves(ss,wfsubTimeRange,sl),'color',colors(ss,:),'linewidth',2.5);
                curBand = find(subPpmAxis>=searchBands(ss,1),1) : ...
                    find(subPpmAxis>=searchBands(ss,2),1);
                plot3(subPpmAxis(curBand),timeAxis(wfsubTimeRange(end))*ones(1,numel(curBand)),...
                    normAvgSpectrum(curBand),'color',colors(ss,:),'linewidth',2.2);
            end
            xlim([ppmRange(1) ppmRange(2)])
            if numel(wfsubTimeRange) > 1
                ylim(timeAxis([wfsubTimeRange(1) wfsubTimeRange(end)]))
            end
            zlim([0 1.1])
            xlabel('Chemical Shift (ppm)');
            ylabel('Time (sec)');
            zlabel('Signal (AU)');
            
            if exist('export_fig', 'file')
                export_fig(curFigWaterfall, fullfile(inDir, sprintf( ...
                    'faFidEcho_waterfallPlot-Slice%d',sl)),'-a1');
            else
                print(curFigWaterfall, fullfile(inDir, sprintf( ...
                    'faFidEcho_waterfallPlot-Slice%d',sl)),'-dpng');
            end
            eval(sprintf('figWaterfall.slice%d = curFigWaterfall;',sl))
            
        end
        
        if numel(timeAxis) > 1
            timeAxis = timeAxis + diff(timeAxis(1:2))/2;
        end
    end
    
end

timeAxis = timeAxis(subTimeRange);  % crop time axis for output

%% Assemble output structure
%{
  spectra     - Complex slice spectra over time
  avgSpectrum - Time-averaged magnitude spectrum
  ppmAxis     - Spectral axis coordinates (ppm)
  timeAxis    - Time axis coordinates (sec)
  timeCurves  - Array of individual species time curves (normalized)
  snrCurves   - Array of individual species SNR curves
  ratioCurves - Array of cumulative ratiometric time curves for  
                each species, first is assumed to be precursor
  chemShiftsPw- Array of point-wise individual species coordinates
                over time (ppm)
  lineWidths  - Vector of individual species linewidths (ppm)
  lineWidthsPw- Array of point-wise individual species linewidths
                over time (ppm)
  spSnr       - SNR of spectrum (combined signal from all species)
  headers     - Combined acqp and method headers
  funVersion  - idFidEcho function version
%}

names = { 'spectra', 'avgSpectrum', 'ppmAxis', 'timeAxis', 'timeCurves', ...
    'snrCurves', 'ratioCurves', 'chemShifts', 'chemShiftsPw', ...
    'lineWidths', 'lineWidthsPw', 'spSnr', 'headers', 'funVersion'};

for nn = 1:numel(names)-2
    eval( sprintf( 'outStruct.%s = %s;', names{nn}, names{nn} ) );
end

if showPlots
    for nn = numel(names)-1:numel(names)
        eval( sprintf( 'outStruct.%s = %s;', names{nn}, names{nn} ) );
    end
end


% save(fullfile(inDir,'autospf_output.mat'),'outStruct');
% keyboard

%% Write nLac values to file and std out

% f = fopen(fullfile(inDir,'autospf_nLac.txt'),'w');
% fprintf(f,fString);
% fclose(f);
fprintf(1, fString);

if isdeployed
    fprintf(1, 'Leaving %s at %s\n', mfilename, datestr(now));
end


return


%####################################################################

function parsein(params)

%{
  searchBands - Array of spectral bands to find peaks (ppm)
      default = [169, 173; 182, 185] (pyr; lac)
  lbHz        - Exponential line broadening factor (Hz)
      default = 10
  qMethod     - String specifying quantification method
      default = 'fwhm'    - FWHM integration
                'genspec' - Fit peak amplitudes in freq domain
  noiseRegion - 2-element vector specifying limits of noise
                region used for baseline correction (ppm)
      default = [195 200]
  timeRegion   - 2-element vector specifying time region for analysis (sec)
      default = [0 150]
  ppmRange    - 2-element vector specifying spectral coordinates for
                cropping waterfall plot (ppm)
      default = [160 190]
  showPlots   - Logical, Show result plots
      default = true
%}

defaults = struct('searchBands', [169 173; 182 185], 'lbHz', 10, ...
    'qMethod', 'fwhm', 'noiseRegion', [195 200], 'timeRegion', [0 150], ...
    'ppmRange', [160 190], 'showPlots', true);

defNames = fieldnames(defaults);
if isstruct(params)
    inNames = fieldnames(params);
    for ii = 1:numel(inNames)
        if ~any( strcmp( inNames{ii}, defNames))
            if ~any( strcmpi( inNames{ii}, defNames))
                warning('faFidEcho:parsein',...
                    'Input parameter "%s" not recognized',inNames{ii})
            else
                ind = find( strcmpi( inNames{ii}, defNames));
                params.(defNames{ind}) = params.(inNames{ii});
            end
        end
    end
end
    
for ii = 1:numel(defNames)
    if ~isfield(params,defNames{ii})
        params.(defNames{ii}) = defaults.(defNames{ii});
    end
    assignin('caller',defNames{ii},params.(defNames{ii}));
end


%####################################################################

function [int,loc,lwPts] = fwhm(spectrum,idxSearch)

spectrum = real(spectrum(:)).';
[pk,loc] = max(spectrum(idxSearch(1):idxSearch(2)));
loc = loc + idxSearch(1) - 1;
locHmLo = loc;  % index of point just above half max in - dir
locHmHi = loc;  % index of point just above half max in + dir

while true  % walk down
    if locHmLo-1 == 1
        fprintf(1, 'FWHM warning: reached first index of spectrum\n');
        break
    elseif spectrum(locHmLo-1) > pk
        break
    elseif spectrum(locHmLo-1) < (pk/2)
        break
    else
        locHmLo = locHmLo - 1;
    end
end

while true  % walk up
    if locHmHi+1 == numel(spectrum)
        fprintf(1, 'FWHM warning: reached last index of spectrum\n');
        break
    elseif spectrum(locHmHi+1) > pk
        break
    elseif spectrum(locHmHi+1) < (pk/2)
        break
    else
        locHmHi = locHmHi + 1;
    end
end

% interpolate to intermediate positions
partLo = interp1( spectrum(locHmLo-[1 0]), [1, 0], pk/2, 'linear');
if isnan(partLo), partLo = 0; end   % low index not at half-max
partHi = interp1( spectrum(locHmHi+[0 1]), [0, 1], pk/2, 'linear');
if isnan(partHi), partHi = 0; end   % high index not at half-max

% fwhm integral
x = locHmLo:locHmHi;
y = spectrum(x);
x = [locHmLo-partLo, x, locHmHi+partHi];
y = [pk/2, y, pk/2];
int = trapz(x,y);
lwPts = x(end) - x(1);

return


%####################################################################

function spectrum = genspec(pars,freqs)

nShifts = numel(pars)/3;
pars = reshape(pars,3,nShifts);
A = pars(1,:);  % amplitude
F = pars(2,:);  % frequency
W = pars(3,:);  % linewidth
N = numel(freqs);

t = (1:N).';
spectrum = zeros(N,1);
for ss = 1:nShifts
    spectrum = spectrum + ( A(ss) * W(ss) ) ./ ( W(ss) + 1i*(t-F(ss)) );
end
spectrum = abs( spectrum );

return
