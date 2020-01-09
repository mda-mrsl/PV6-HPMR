function [g, traj] = create2dEpi(epiType, fov, mtx, dw, Gmax, Smax, tGrad, gmr, pfy, saveFile)
%CREATE2DEPI designs 2D flyback or symmetric EPI trajectories for snapshot imaging
%
%   Usage: [g, traj] = create2dEpi(epiType, fov, mtx, dw, Gmax, Smax, tGrad, gmr, pfy, saveFile)
%          [g, traj] = create2dEpi(optStruct)
%
%   Inputs:
%     epiType,  false for flyback, true for symmetric (default false)
%     fov,      Field of view (default 4)               [cm]
%     mtx,      Matrix size (default 16)
%     dw,       Acquisition dwell time (default 32)     [us]
%     Gmax,     Max gradient amplitude (default 600)    [mT/m]
%     Smax,     Max gradient slew rate (default 4e3)    [T/m/s]
%     tGrad,    Gradient timing resolution (default 8)  [us]
%     gmr,      Gyromagnetic ratio (default 10.71, 13C) [MHz/T]
%     pfy,      Fraction of k-space sampled in blipped direction (default 1)
%     saveFile, String specifying file to write acquisition time values to 
%               as text for recon, if false (default) no file is written.
%               If -1, output files similar to create2dEpi_codegen are saved.
% 
%     optStruct, Structure containing above inputs as field/value pairs.
%                Default values are used for omitted inputs.
%                Additional fields may be specified for this input syntax:
%       showPlot, If true plot trajectory (default true)
%       t2,       T2* relaxation time constant for PSF (default 20) [ms]
%       showPsf,  If true display point spread function (default false)
%                
%   Outputs:
%     g,    2D snapshot EPI readout gradient in form Gx+iGy    [mT/m]
%     traj, trajectory structure containing fields:
%       t,      Gradient time vector                             [ms]
%       g,      2D snapshot EPI readout gradient in form Gx+iGy  [mT/m]
%       T,      Acquisition time vector (values >0 are samples)  [ms]
%       K,      Acquisition K-space points                       [cycles/cm]
%       idxAcq, Acquisition index mask
%       TE,     Effective echo time                              [ms]
%       ESP,    Echo spacing                                     [ms]
%       encEff, Encoding efficiency (fraction of time acquiring samples)
%       matrix, Acquired matrix size (mtx x mty)
%       inputs, Inputs passed to function (optStruct filled with defaults)
%  
%   See also FIND2DEPIPARAMS, SHAPEBRUKERFIDALL, RECOBRUKERKSPACE
% 
% 06/2019, Keith Michel

%% Parse Inputs
scriptVersion = 'v20190923';

if ~nargin, help(mfilename); return; end

if isstruct(epiType)
    optStruct = epiType;
    parseIn(optStruct);
else
    if nargin < 10, help(mfilename); return; end
    t2       = false;
    showPlot = false;
    showPsf  = false;
end

fprintf('Entering %s (%s)\n\n', mfilename, scriptVersion);

tol = 1e-3; % roundoff tolerance for acquisition window k-vals
g   = 0+0i;
k   = 0+0i;
epiType = epiType == 1;

if rem(dw, tGrad*2)
    error('create2dEpi:dwellTime', ...
        'Acquisition dwell time must be even integer multiple of gradient timing resolution');
else
    np = round(dw/tGrad);
    dw = tGrad * np;
end

mtx = double(mtx);
pfy = max(min(pfy, 1), 0.5);
mty = ceil(mtx*pfy);
pfy = mty/mtx;

optStruct = struct('epiType', epiType, 'fov', fov, 'mtx', mtx, ...
    'dw', dw, 'Gmax', Gmax, 'Smax', Smax, 'tGrad', tGrad, 'gmr', gmr, ...
    'pfy', pfy, 'saveFile', saveFile, 't2', t2, 'showPlot', showPlot, ...
    'showPsf', showPsf);

if ~epiType
    fprintf('Flyback EPI, %g cm FOV, %d x %d mtx, %g Hz BW \n', fov, int16(mtx), int16(mty), 1e6/dw);
else
    fprintf('Symmetric EPI, %g cm FOV, %d x %d mtx, %g Hz BW \n', fov, int16(mtx), int16(mty), 1e6/dw);
end

%% Design readout gradient
% force ramp duration to be integer multiple of dw/2
% for symmetric waveforms, force phase blip duration to be same as ramps
dk    = 1/fov;                                  % [cycles/cm]
gRead = 1e5 * dk / gmr / dw;                    % [mT/m]
assert(gRead <= Gmax, 'create2dEpi:maxGrad', ...
        'Maximum gradient amplitude exceeded on readout. Increase FOV or dwell time.');

nIter = 0;
fRamp = 0.5;
gRamp = 0;
if ~epiType
    fprintf('Designing readout ramp');
else
    fprintf('Designing readout ramp and phase blip');
    gBlip = 0;
end
while true
    fprintf('.');
    nIter = nIter + 1;
    assert(nIter < 40, 'create2dEpi:rampIter', ...
        'Maximum iterations met in ramp design while loop');
    nRamp = fRamp*np;
    gTmp  = (0:nRamp)*gRead/nRamp;
    sRamp = diff(gTmp(1:2))*1e3/tGrad;
    if ~epiType
        if sRamp > Smax
            fRamp = fRamp + 0.5;
        else
            gRamp = gTmp;
            break;
        end
    else
        gBlip = [gTmp, fliplr(gTmp)];
        kBlip = trapz(gBlip) * gmr * tGrad * 1e-5;
        gBlip = gBlip * dk / kBlip;
        sBlip = diff(gBlip(1:2)) * 1e3 / tGrad;
        if sRamp > Smax || sBlip > Smax || max(gBlip) > Gmax
            fRamp = fRamp + 0.5;
        else
            gRamp = gTmp;
            break;
        end
    end
end
fprintf('\n\t%g us readout ramp, %g acquisition points\n', nRamp*tGrad, fRamp);
kRamp = trapz(gRamp) * gmr * tGrad * 1e-5;                  % [cycles/cm]
tRead = dw * mtx;                                           % [us]
nRead = round(tRead/tGrad);
gRead = [gRamp, gRead*ones(1,nRead), fliplr(gRamp)];
kRead = trapz(gRead) * gmr * tGrad * 1e-5;                  % [cycles/cm]
fprintf('\t%g us readout, %g mT/m gradient\n', tRead, max(gRead));

%% Design flyback/prephasing gradient
% force duration to be integer multiple of dw
nIter = 0;
fRise = 0.5;
fPlat = 1;
kMax  = mtx / fov / 2;
if ~epiType
    fprintf('Designing flyback gradient');
    gFlyb = 0;
else
    fprintf('Designing prephasing gradient');
    gPre = 0;
end
while true
    fprintf('.');
    nIter = nIter + 1;
    assert(nIter < 40, 'create2dEpi:flybIter', ...
            'Maximum iterations met in flyback lobe design while loop');
    nRise = fRise*np;
    nPlat = fPlat*np;
    gTmp  = (0:nRise)*Smax*tGrad/1e3;
    gTmp(gTmp>Gmax) = Gmax;
    gTmp  = [gTmp, gTmp(end)*ones(1,nPlat), fliplr(gTmp)];
    kTmp  = trapz(gTmp) * gmr * tGrad * 1e-5;
    if ~epiType
        if kTmp >= kRead
            gFlyb = gTmp * kRead / kTmp;
            break;
        end
    else
        if kTmp >= (kMax+kRamp+dk)
            gPre = gTmp * (kMax+kRamp+dk) / kTmp;
            break;
        end
    end
    if max(gTmp) < Gmax
        fRise = fRise + 0.5;
    end
    fPlat = fPlat + 1;
end
if ~epiType
    fprintf('\n\t%g us flyback gradient, %g acquisition points\n', ...
        numel(gFlyb)*tGrad, 2*fRise+fPlat);
else
    fprintf('\n\t%g us prephasing gradient, %g acquisition points\n', ...
        numel(gPre)*tGrad, 2*fRise+fPlat);
end

%% Design/validate prephasing gradients
% for flyback, scale flyback grad
% for symmetric, scale prephasing grad
if ~epiType
    gPreY   = gFlyb * (kMax - dk*(mtx-mty)) / kRead;
    gPreX   = gFlyb * (kMax+kRamp+dk) / kRead;
    gBlip   = gFlyb * dk / kRead;
    nEchoes = mty;
else
    gPreY   = gPre * (kMax - dk*(mtx-mty)) / (kMax+kRamp+dk);
    gPreX   = gPre;
    nEchoes = ceil(mty/2);
end

% Adjust read prephasing gradient to correct acquisition k-values
nzPre = np;   % 1 acq pt, is there a better value to use here?
gTmp  = [-gPreX, zeros(1,nzPre), gRead(2:end)];
t     = tGrad *(0:numel(gTmp)-1);
T     = dw*(0:floor(t(end)/dw));
kTmp  = cumtrapz(gTmp) * gmr * tGrad * 1e-5;
% It works best to use 1/10 tol here for rounding KTmp. Not sure why...
KTmp  = round(interp1(t,kTmp,T)/tol*10)*tol/10;
idx   = find(fliplr(KTmp) <= kMax-dk, 1) - 1;
kErr  = KTmp(end-idx(1)) - kMax+dk;
gPreX = gPreX * (kMax+kRamp+dk+kErr) / (kMax+kRamp+dk);

%% Assemble gradient waveform
% insert zeros between flybacks and readouts until k-values are correct
nz     = 0;
kFind  = -kMax:dk:kMax-dk;
idxADC = 0;
fprintf('Assembling full gradient waveform');
while true
    fprintf('.');
    assert(nz < 1+dw/tGrad, 'create2dEpi:idxAcq', ...
            '\nFailed to assemble gradient waveform. Can''t find acquisition window.');
    if ~epiType
        echoTrain = repmat([gRead(2:end), -gFlyb(2:end), zeros(1,nz)], 1, mty-1) + ...
            1i*repmat([zeros(1,numel(gRead)-1), gBlip(2:end), zeros(1,nz)], 1, mty-1);
        g = [-gPreX-1i*gPreY, zeros(1,nzPre), echoTrain, gRead(2:end)];
    else
        nEchoes    = ceil(mty/2);
        echoTrainX = repmat([gRead(2:end), -gRead, zeros(1,nz)], 1, floor(mty/2));
        if rem(mty,2)
            echoTrainX = [echoTrainX, gRead(2:end)];
            echoTrainY = [zeros(1,nRamp+nRead), gBlip, ...
                repmat([zeros(1,nRead+nz), gBlip, zeros(1,nRead-1), gBlip], 1, floor((mty-1)/2)-1)];
            echoTrainY = [echoTrainY, zeros(1,nRead-1), gBlip];
        else
            echoTrainY = [zeros(1,nRamp+nRead), gBlip, ...
                repmat([zeros(1,nRead+nz), gBlip, zeros(1,nRead-1), gBlip], 1, floor((mty-1)/2))];
        end
        echoTrain = echoTrainX + ...
            1i*[echoTrainY, zeros(1,numel(echoTrainX)-numel(echoTrainY))];
        g = [-gPreX-1i*gPreY, zeros(1,nzPre), echoTrain];
    end
    t = 1e-3*tGrad*(0:numel(g)-1);
    T = 1e-3*dw*(0:floor(1e3*t(end)/dw));
    k = (cumtrapz(real(g)) + 1i*cumtrapz(imag(g))) * gmr * tGrad * 1e-5;
    K = round(interp1(t, real(k), T) / tol) * tol ;
    idx = findPattern2(K, round(kFind / tol) * tol);
    
    if nnz(idx) < nEchoes
        nz = nz+1;
    else
        idxADC = idx;
        break;
    end
end

if ~epiType
    g = [g, -gFlyb*real(k(end))/kRead - 1i*gFlyb*imag(k(end))/kRead];
else
    g = [g, -gPre*real(k(end))/(kMax+kRamp+dk) - 1i*gPre*imag(k(end))/(kMax+kRamp+dk)];
end
t = 1e-3*tGrad*(0:numel(g)-1);
T = 1e-3*dw*(0:floor(1e3*t(end)/dw));
k = (cumtrapz(real(g)) + 1i*cumtrapz(imag(g))) * gmr * tGrad * 1e-5;
K = interp1(t,real(k),T) + 1i*interp1(t,imag(k),T);
idxAcq  = false(size(T));
idxDiff = unique(diff(idxADC(idxADC>0)));
assert(numel(idxDiff) == 1, 'create2dEpi:idxDiff', ...
    'Echo spacing is not uniform in time');
idxDiff = idxDiff(1);
for ii=1:nEchoes
    idxAcq(idxADC(ii):idxADC(ii)+mtx-1) = true;
    if epiType
        idxAcq(idxADC(ii)+ceil(idxDiff/2):idxADC(ii)+ceil(idxDiff/2)+mtx-1) = true;
    end
end
if epiType && rem(mty,2)
    idxAcq(idxADC(ii)+ceil(idxDiff/2):idxADC(ii)+ceil(idxDiff/2)+mtx-1) = false;
    idxAcq(numel(T)+1:end) = [];
end
tAcq    = T(idxAcq);
kAcq    = interp1(t,real(k),tAcq) + 1i*interp1(t,imag(k),tAcq);
[~,idx] = min(abs(kAcq));
TE      = tAcq(idx(1));
if ~epiType
    ESP = tAcq(mtx+1)-tAcq(1);
else
    [~,idx] = min(abs(kAcq(1:mtx)));
    echo1   = tAcq(idx(1));
    [~,idx] = min(abs(kAcq(mtx+1:mtx*2)));
    echo2   = tAcq(idx(1)+mtx);
    ESP     = echo2 - echo1;
end
encEff  = numel(tAcq)/numel(T);
fprintf(['\nTotal duration %g ms, TE %g ms\n', ...
    '\tESP %g ms, Phase Bandwidth %g Hz\n\tEncoding Efficiency %g%%\n'], ...
    t(end), TE, ESP, 1e3/ESP, 100*encEff);

% Acquisition time values > 0 are recon samples
T(~idxAcq) = -T(~idxAcq);


%% Assemble traj output structure
traj = struct('t', t, 'g', g, 'T', T, 'K', K, 'idxAcq', idxAcq, ...
    'TE', TE, 'ESP', ESP, 'encEff', encEff, 'matrix', [mtx, mty], ...
    'inputs', optStruct);

%% Save acquisition time values to text for recon
if ~saveFile
    fprintf('Not saving acquisition times to file\n');
else
    if ~ischar(saveFile)
        gmrName = num2str(round(gmr*100)/100);
        if strcmpi(gmrName, '10.71')
            gmrName = '13C';
        elseif strcmpi(gmrName, '42.58')
            gmrName = '1H';
        else
            gmrName = 'XX';
        end
        if ~epiType
            epiName = 'fbEPI';
        else
            epiName = 'symEPI';
        end
        dwRounded = '';
        if dw ~= round(dw)
            dwRounded = 'R';
        end
        if saveFile ~= -1
            saveFile = sprintf('%s_%s_fov%dmm_%dx%d_dw%dus%s_%s_%s', ...
                epiName, gmrName, int16(fov*10), int16(mtx), int16(mty), ...
                int16(dw), dwRounded, mfilename, scriptVersion);
        end
    end
    if saveFile == -1
        fileName1 = sprintf('tAcq_%s_%s_fov%dmm_%dx%d_dw%dus%s_%s_%s', ...
            epiName, gmrName, int16(fov*10), int16(mtx), int16(mty), ...
            int16(dw), dwRounded, mfilename, scriptVersion);
        fileName2 = sprintf('gradWave_%s_%s_fov%dmm_%dx%d_dw%dus%s_%s_%s.gp', ...
            epiName, gmrName, int16(fov*10), int16(mtx), int16(mty), ...
            int16(dw), dwRounded, mfilename, scriptVersion);
        fprintf(1, 'Saving acquisition times to %s and gradient waveform to %s\n', ...
            fileName1, fileName2);
        fid1 = fopen(fileName1, 'w');
        fprintf(fid1, ['# Acquisition timepoints in milliseconds from %s (%s). ', ...
            'Positive values indicate recon samples\n'], mfilename, scriptVersion);
        fprintf(fid1, ['# Inputs: FOV = %g [cm], MTX = %d, DW = %g [us], ', ...
            'GMAX = %g [mT/m], SMAX = %g [T/m/s], TGRAD = %g [us], ', ...
            'GMR = %g [MHz/T], PFY = %g, SAVEFILE = %s\n'], ...
            fov, int16(mtx), dw, Gmax, Smax, tGrad, gmr, pfy, saveFile);
        fprintf(fid1, ['# Waveform Parameters: TOTALTIME = %g [ms], ', ...
            'ECHOTIME = %g [ms], ECHOSPACING = %g [ms], ', ...
            'PHASEBW = %g [Hz], ENCODINGEFFICIENCY = %g %%\n'], ...
            t(end), TE, ESP, 1e3/ESP, 100*encEff);
        for ii = 1:numel(T)
            fprintf(fid1, '%.6f \n', T(ii));
        end
        fid2 = fopen(fileName2, 'w');
        for ii = 1:numel(t)
            fprintf(fid2, '%.6f %.6f %.6f\n', real(g(ii)), imag(g(ii)), 0);
        end
    else
        fprintf('Saving gradient waveform in [mT/m] to %s \n', saveFile);
        fid = fopen([saveFile, '.gp'], 'w');
        for ii = 1:numel(g)
            fprintf(fid, '%.6e %.6e %6e\n', real(g(ii)), imag(g(ii)), 0);
        end
        save(saveFile, 'g', 'traj', 'optStruct');
    end
    fclose('all');
end

%% Plot trajectory
if showPlot
    figure('position', [100 100 1200 800]);
    colors = lines(2);
    
    subplot(3, 2, 1)
    h1 = area(t, real(g), 'EdgeColor', colors(1,:), 'FaceColor', colors(1,:), ...
        'FaceAlpha', 0.2); 
    hold on
    h2 = plot(tAcq, interp1(t,real(g),tAcq), 'o', 'MarkerEdgeColor', colors(2,:));
    xlabel('Time (ms)'); 
    ylabel('Gx (mT/m)');
    set(gca, 'fontsize', 15); 
    grid on
    
    subplot(3, 2, 2)
    h1 =area(t, imag(g), 'EdgeColor', colors(1,:), 'FaceColor', colors(1,:), ...
        'FaceAlpha', 0.2); 
    hold on
    h2 = plot(tAcq, interp1(t,imag(g),tAcq), 'o', 'MarkerEdgeColor', colors(2,:));
    xlabel('Time (ms)'); 
    ylabel('Gy (mT/m)');
    set(gca, 'fontsize', 15); 
    grid on
    
    subplot(3, 2, [3, 5])
    h1 = plot(real(k), imag(k), '.-'); 
    hold on
    h2 = plot(real(kAcq), imag(kAcq), 'o');
    xlabel('Kx (cycles/cm)'); 
    ylabel('Ky (cycles/cm)');
    legend([h1, h2], 'Gradient', 'Samples', 'location', 'northoutside', ...
        'orientation', 'horizontal')
    set(gca, 'fontsize', 15); 
    grid on
    axis(1.1*kMax * [-1 1 -1 1])
    axis ij square
    
    subplot(3, 2, [4,6])
    h1 = plot3(t, real(k), imag(k), '.-'); 
    hold on
    h2 = plot3(tAcq, real(kAcq), imag(kAcq), 'o');
    xlabel('Time (ms)'); 
    ylabel('Kx (cycles/cm)'); 
    zlabel('Ky (cycles/cm)');
    ylim(1.1*kMax * [-1 1]); 
    zlim(1.1*kMax * [-1 1]);
    set(gca, 'fontsize', 15); 
    grid on
    axis ij 
    
end

%% Plot point spread function
if showPsf  
    sig = exp(-tAcq/t2);
    psf = ifftdim(reshape(sig, [mtx, mty]), 1:2).';
    figure;
    posx = linspace(-fov*5, fov*5, mtx);
    posy = linspace(-fov*5, fov*5, mty);
    imagesc(posx, posy, abs(psf));
    xlabel('Readout (mm)')
    ylabel('Blipped (mm)')
    set(gca, 'fontsize', 15)
    axis square
    title({'Point spread function', sprintf('T2* = %g ms', t2)}, 'fontsize', 16)
    if ~verLessThan('matlab', '9.3')
        set(gca, 'colorscale', 'log')
        cb = colorbar();
        cb.Ruler.Scale = 'log';
        cb.Ruler.MinorTick = 'on';
    end
end


fprintf('\nLeaving %s (%s)\n', mfilename, scriptVersion);

end

%####################################################################

function parseIn(opts)

defaults = struct('epiType', false, 'fov', 4, 'mtx', 16, 'dw', 32, ...
    'Gmax', 600, 'Smax', 4e3, 'tGrad', 8, 'gmr', 10.71, 'pfy', 1, ...
    't2', 20, 'saveFile', false, 'showPlot', true, 'showPsf', false);

names = fieldnames(defaults);
for ii = 1:numel(names)
    if ~isfield(opts, names{ii})
            opts.(names{ii}) = defaults.(names{ii});
    end
    assignin('caller', names{ii},opts.(names{ii}));
end

end

function start = findPattern2(array, pattern)
%findPattern2 Locate a pattern in an array.
%
%   indices = findPattern2(array, pattern) finds the starting indices of
%   pattern within array.
%
%   Example:
%   a = [0 1 4 9 16 4 9];
%   patt = [4 9];
%   indices = findPattern2(a,patt)
%   indices =
%        3     6

% Let's assume for now that both the pattern and the array are non-empty
% VECTORS, but there's no checking for this. 
% For this algorithm, I loop over the pattern elements.
len = length(pattern);
% First, find candidate locations; i.e., match the first element in the
% pattern.
start = find(array==pattern(1));
% Next remove start values that are too close to the end to possibly match
% the pattern.
endVals = start+len-1;
start(endVals>length(array)) = [];
% Next, loop over elements of pattern, usually much shorter than length of
% array, to check which possible locations are valid still.
for pattval = 2:len
    % check viable locations in array
    locs = pattern(pattval) == array(start+pattval-1);
    % delete false ones from indices
    start(~locs) = [];
end

end
