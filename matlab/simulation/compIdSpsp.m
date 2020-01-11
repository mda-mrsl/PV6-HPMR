% This script simulates relative SNR and NRMSE for spectral encoding via
% spectral-spatial excitation and IDEAL encoding (data for Figure 5).
% 
% Please note, as written these simulations will take a considerable amount
% of time.

clear variables

trajFile = 'flybackEPI_fov40_mtx16_dw48_13C_20190103';
load(['../../hpMR/gp/', trajFile], 'mtx', 'tAcq')

rng('default')

%% Create phantom, hard-coded for mtx=16
% [val,	 xwd,	 ywd,   xc,	   yc,	ang]
E = [1,	0.25,   0.25, -0.4,	 0.65, 0; ...
     2,	0.25,   0.25, -0.4,     0, 0; ...
     3,	0.25,   0.25, -0.4, -0.65, 0; ...
     4,	0.25,   0.25,  0.4,  0.65, 0; ...
     5,	0.25,   0.25,  0.4,     0, 0; ...
     6,	0.25,   0.25,  0.4, -0.65, 0];
nMet = 6
P = phantom(mtx, E);

%% Initialize variables
% Test variables
fa     = [2:2:90];       % [deg], effective for Mz loss (IDEAL fa are less)
t2s    = [2:2:100];      % [ms]
noiVar = [0.01:0.01:0.5]; % max image SNR equals mtx/sqrt(noiVar) (see snrCalc.m)
% fa     = [2:8:88];  % [deg], effective for Mz loss (IDEAL fa are less)
% t2s    = [2:8:50];  % [ms]
% noiVar = [0.01:.05:0.5];

% fa     = [10 90];  % [deg], effective for Mz loss
% t2s    = [20 1e9];
% noiVar = [0.1 1];

% nExc  = 6:16;      % # of IDEAL excitations
% dte   = [0.46, 0.48, 0.53, 0.48, 0.43, 0.40, 0.37, 0.49, 0.45, 0.42, 0.53] ;   % [ms]
nExc  = 8;      % # of IDEAL excitations
dte   = 0.52;   % [ms]
nNexc = numel(nExc);

% Indices of images to keep
imFa  = 5:5:45;       % 10:10:90 deg
imT2s = 5:5:50;       % 10:10:100 ms
imNoi = [1, 5:5:50];  % 0.01, 0.05:0.05:0.5
imRep = 5;            % number of imgs per iter to keep
% imFa  = 2:2:10;
% imT2s = 2:2:6;
% imNoi = 2:2:10;

% Intermediate variables
nReps = 100
% nReps = 10
fHz   = [921, 642, 453, 0, -566, -755]; % LHAPUB
tSp   = tAcq - tAcq(mtx^2/2 + mtx/2 + 1) + 13.990;  % [ms]
tId   = tAcq - tAcq(mtx^2/2 + mtx/2 + 1) + 11.502;  % [ms]
phT1  = 700; % [ms]
hpT1  = 1e4; % [ms]
idTr  = 50;  % [ms]

% Output variables
resSize  = [numel(fa), numel(t2s), numel(noiVar)];
nIter    = prod(resSize)
imSize   = [numel(imFa), numel(imT2s), numel(imNoi)];
imIter   = prod(imSize)

spSnrOutMean = zeros(nIter, nMet);
spSnrOutStd  = zeros(nIter, nMet);
spSnrOutMin  = zeros(nIter, nMet);
spSnrOutMax  = zeros(nIter, nMet);

idSnrOutMean = zeros(nIter, nNexc, nMet);
idSnrOutStd  = zeros(nIter, nNexc, nMet);
idSnrOutMin  = zeros(nIter, nNexc, nMet);
idSnrOutMax  = zeros(nIter, nNexc, nMet);

idPhSnrOutMean = zeros(nIter, nNexc, nMet);
idPhSnrOutStd  = zeros(nIter, nNexc, nMet);
idPhSnrOutMin  = zeros(nIter, nNexc, nMet);
idPhSnrOutMax  = zeros(nIter, nNexc, nMet);

spNrmseMean = zeros(nIter, nMet);
spNrmseStd  = zeros(nIter, nMet);
spNrmseMin  = zeros(nIter, nMet);
spNrmseMax  = zeros(nIter, nMet);

idNrmseMean = zeros(nIter, nNexc, nMet);
idNrmseStd  = zeros(nIter, nNexc, nMet);
idNrmseMin  = zeros(nIter, nNexc, nMet);
idNrmseMax  = zeros(nIter, nNexc, nMet);

idPhNrmseMean = zeros(nIter, nNexc, nMet);
idPhNrmseStd  = zeros(nIter, nNexc, nMet);
idPhNrmseMin  = zeros(nIter, nNexc, nMet);
idPhNrmseMax  = zeros(nIter, nNexc, nMet);

spNrmseNf    = zeros(nIter, nMet);
idNrmseNf    = zeros(nIter, nNexc, nMet);
idPhNrmseNf  = zeros(nIter, nNexc, nMet);

spIms     = zeros(mtx, mtx, nIter, imRep, nMet);
idIms     = zeros(mtx, mtx, nIter, imRep, nNexc, nMet);
idPhIms   = zeros(mtx, mtx, nIter, imRep, nNexc, nMet);
spImsNf   = zeros(mtx, mtx, nIter, nMet);
idImsNf   = zeros(mtx, mtx, nIter, nNexc, nMet);
idPhImsNf = zeros(mtx, mtx, nIter, nNexc, nMet);

%% Iterate through
rng('default')
input('Press Enter to continue.');
startTime = datestr(now, 30);
fprintf('Starting %d iterations at %s.\n', nIter, startTime)
tic
flt = ones(mtx);
% flt = hann(mtx);
% flt = (flt .* flt.').^2;
parfor aa = 1:nIter
    [a,b,c] = ind2sub(resSize, aa);
    curFa   = fa(a);
    curT2s  = t2s(b);
    curNvar = noiVar(c);
    useImgs = all([any(imFa == a), any(imT2s == b), any(imNoi == c)]);
    curSpSnrOut    = zeros(nReps, nMet);
    curIdSnrOut    = zeros(nReps, nNexc, nMet);
    curIdPhSnrOut  = zeros(nReps, nNexc, nMet);
    curSpNrmse     = zeros(nReps, nMet);
    curIdNrmse     = zeros(nReps, nNexc, nMet);
    curIdPhNrmse   = zeros(nReps, nNexc, nMet);
    curSpNrmseNf   = zeros(nMet, 1);
    curIdNrmseNf   = zeros(nNexc, nMet);
    curIdPhNrmseNf = zeros(nNexc, nMet);
    if useImgs
        curSpIms       = zeros(mtx, mtx, imRep, nMet);
        curIdIms       = zeros(mtx, mtx, imRep, nNexc, nMet);
        curIdPhIms     = zeros(mtx, mtx, imRep, nNexc, nMet);
        curSpImsNf     = zeros(mtx, mtx, nMet);
        curIdImsNf     = zeros(mtx, mtx, nNexc, nMet);
        curIdPhImsNf   = zeros(mtx, mtx, nNexc, nMet);
    end
    
    tmpId   = cell(1, nNexc);
    tmpIdPh = cell(1, nNexc);
    for nn = 1:nNexc
        tmpId{nn}   = zeros(mtx, mtx, nExc(nn));
        tmpIdPh{nn} = zeros(mtx, mtx, nExc(nn));
    end
    
    for bb = 1:nMet
        sigMask = P==bb;
        k = flt .* fftdim(sigMask);
        p = double(sigMask);
        
        % Spectral-spatial excitation w perfect excitation profile:
        kSp = k * sind(curFa) .* reshape(exp(-tSp/curT2s), mtx, mtx);
        mSp = abs(ifftdim(kSp));
        
        kNoise = sqrt(curNvar) * (randn(mtx, mtx, nReps) + 1i*randn(mtx, mtx, nReps));
        mNoise = abs(ifftdim(kNoise, 1:2));
        kSpNsy = repmat(kSp, 1, 1, nReps) + kNoise;
        mSpNsy = abs(ifftdim(kSpNsy, 1:2));
        tmp      = nan(mtx, mtx, nReps);
        msk      = repmat(sigMask, 1, 1, nReps);
        tmp(msk) = mSpNsy(msk);
        tmp      = reshape(tmp, mtx^2, nReps);
        mNoise   = reshape(mNoise, mtx^2, nReps);
        curSpSnrOut(:,bb) = sqrt(2-pi/2) * mean(tmp, 1, 'omitnan') ./ std(mNoise, 0, 1);
        
        for ii = 1:nReps
            tmp = nscale(mSpNsy(:,:,ii));
            curSpNrmse(ii,bb) = norm(tmp(:) - p(:))/norm(p(:));
        end
        
        if useImgs, curSpIms(:,:,:,bb) = mSpNsy(:,:,1:imRep); end
        if c == 1
            if useImgs, curSpImsNf(:,:,bb)  = mSp; end
            curSpNrmseNf(bb) = norm(nscale(mSp(:)) - p(:))/norm(p(:));
        end
        
        % IDEAL
        for nn = 1:nNexc
            curNExc = nExc(nn);
            curFaId = vfaConstant(curNExc, curFa);
            kIdHp   = zeros(mtx, mtx, curNExc);
            kIdPh   = zeros(mtx, mtx, curNExc);
            mzHp    = 1;
            mzTh    = 1;
            for cc = 1:curNExc
                curTe = tId+(cc-1)*dte(nn);
                
                % HP system, signal depletes with T1
                kIdHp(:,:,cc) = mzHp * sind(curFaId(cc)) * k .* ...  % exc angle
                    reshape(exp(1i*2*pi*fHz(bb)*curTe/1e3) .* ... % chem shift
                    exp(-curTe/curT2s), mtx, mtx);                % t2*
                mzHp  = mzHp * cosd(curFaId(cc)) * exp(-idTr/hpT1);
                
                % Thermal system, signal recovers with T1
                kIdPh(:,:,cc) = mzTh * sind(curFaId(cc)) * k .* ...  % exc angle
                    reshape(exp(1i*2*pi*fHz(bb)*curTe/1e3) .* ... % chem shift
                    exp(-curTe/curT2s), mtx, mtx);                % t2*
                mzTh  = 1 - (1-mzTh*cosd(curFaId(cc))) * exp(-idTr/phT1);
            end
            tmpId{nn}   = tmpId{nn} + kIdHp;
            tmpIdPh{nn} = tmpIdPh{nn} + kIdPh;
        end
    end
    spSnrOutMean(aa,:) = reshape(mean(curSpSnrOut, 1), [1,nMet]);
    spSnrOutStd(aa,:)  = reshape(std(curSpSnrOut, [], 1), [1,nMet]);
    spSnrOutMin(aa,:)  = reshape(min(curSpSnrOut, [], 1), [1,nMet]);
    spSnrOutMax(aa,:)  = reshape(max(curSpSnrOut, [], 1), [1,nMet]);
    spNrmseMean(aa,:)  = reshape(mean(curSpNrmse, 1), [1,nMet]);
    spNrmseStd(aa,:)   = reshape(std(curSpNrmse, [], 1), [1,nMet]);
    spNrmseMin(aa,:)   = reshape(min(curSpNrmse, [], 1), [1,nMet]);
    spNrmseMax(aa,:)   = reshape(max(curSpNrmse, [], 1), [1,nMet]);
    spNrmseNf(aa,:)    = reshape(curSpNrmseNf, 1, nMet);
    if useImgs
        spIms(:,:,aa,:,:)   = reshape(curSpIms, mtx, mtx, 1, imRep, nMet);
        spImsNf(:,:,aa,:) = reshape(curSpImsNf, mtx, mtx, 1, nMet);
    end
    
    for nn = 1:nNexc
        curNExc = nExc(nn);
        curDte  = dte(nn);
        ideal = @(dat) idDecomp(reshape(dat, mtx^2, curNExc), ...
            fHz, (0:curNExc-1)*curDte, tAcq, curT2s);
        
        kIdHp   = tmpId{nn};
        kDcmpNf = ideal(reshape(kIdHp, mtx^2, curNExc));
        kDcmpNf = reshape(kDcmpNf, mtx, mtx, nMet);

        kIdPh   = tmpIdPh{nn};
        kPhDcmpNf = ideal(reshape(kIdPh, mtx^2, curNExc));
        kPhDcmpNf = reshape(kPhDcmpNf, mtx, mtx, nMet);
        
        for bb = 1:nMet
            sigMask = P==bb;
            p = double(sigMask);
            
            mId   = abs(ifftdim(kDcmpNf(:,:,bb)));            
            mIdPh = abs(ifftdim(kPhDcmpNf(:,:,bb)));
            if c == 1
                if useImgs, curIdImsNf(:,:,nn,bb) = mId; end
                curIdNrmseNf(nn,bb)     = norm(nscale(mId(:)) - p(:))/norm(p(:));
                if useImgs, curIdPhImsNf(:,:,nn,bb) = mIdPh; end
                curIdPhNrmseNf(nn,bb)   = norm(nscale(mIdPh(:)) - p(:))/norm(p(:));
            end
        end
        
        kNoise = sqrt(curNvar) * (randn(mtx, mtx, curNExc, nReps) + ...
            1i*randn(mtx, mtx, curNExc, nReps));
        kIdNsy = repmat(kIdHp, 1, 1, 1, nReps) + kNoise;
        kPhNoise = sum(sqrt(curNvar) * (randn(mtx,mtx,curNExc,nReps,nMet) + ...
                1i*randn(mtx,mtx,curNExc,nReps,nMet)), 5);
        kIdPhNsy = repmat(kIdPh, 1, 1, 1, nReps) + kPhNoise;
        kNDcmp   = zeros(mtx, mtx, nMet, nReps);
        kDcmp    = zeros(mtx, mtx, nMet, nReps);
        kPhNDcmp = zeros(mtx, mtx, nMet, nReps);
        kPhDcmp  = zeros(mtx, mtx, nMet, nReps);
        for ii = 1:nReps
            kNDcmp(:,:,:,ii)   = reshape(ideal(kNoise(:,:,:,ii)), mtx, mtx, nMet);
            kDcmp(:,:,:,ii)    = reshape(ideal(kIdNsy(:,:,:,ii)), mtx, mtx, nMet);
            kPhNDcmp(:,:,:,ii) = reshape(ideal(kPhNoise(:,:,:,ii)), mtx, mtx, nMet);
            kPhDcmp(:,:,:,ii)  = reshape(ideal(kIdPhNsy(:,:,:,ii)), mtx, mtx, nMet);
        end
        mNoise   = abs(ifftdim(kNDcmp, 1:2));
        mIdNsy   = abs(ifftdim(kDcmp, 1:2));
        mPhNoise = abs(ifftdim(kPhNDcmp, 1:2));
        mIdPhNsy = abs(ifftdim(kPhDcmp, 1:2));

        for bb = 1:nMet
            sigMask = P==bb;
            p = double(sigMask);
            
            msk = repmat(sigMask, 1, 1, nReps);
            tmp = nan(mtx, mtx, nReps);
            tms = squeeze(mIdNsy(:,:,bb,:));
            tmn = squeeze(mNoise(:,:,bb,:));
            tmp(msk) = tms(msk);
            tmp = reshape(tmp, mtx^2, nReps);
            tmn = reshape(tmn, mtx^2, nReps);
            curIdSnrOut(:,nn,bb) = sqrt(2-pi/2) * ...
                mean(tmp, 1, 'omitnan') ./ std(tmn, 0, 1);
            if useImgs, curIdIms(:,:,:,nn,bb) = tms(:,:,1:imRep); end
            tmp = nan(mtx, mtx, nReps);
            tms = squeeze(mIdPhNsy(:,:,bb,:));
            tmn = squeeze(mPhNoise(:,:,bb,:));
            tmp(msk) = tms(msk);
            tmp = reshape(tmp, mtx^2, nReps);
            tmn = reshape(tmn, mtx^2, nReps);
            curIdPhSnrOut(:,nn,bb) = sqrt(2-pi/2) * ...
                mean(tmp, 1, 'omitnan') ./ std(tmn, 0, 1);
            if useImgs, curIdPhIms(:,:,:,nn,bb) = tms(:,:,1:imRep); end
            
            for ii = 1:nReps
                tmp = nscale(mIdNsy(:,:,bb,ii));
                curIdNrmse(ii,nn,bb) = norm(tmp(:) - p(:))/norm(p(:));
                tmp = nscale(mIdPhNsy(:,:,bb,ii));
                curIdPhNrmse(ii,nn,bb) = norm(tmp(:) - p(:))/norm(p(:));
            end

        end

    end
    idSnrOutMean(aa,:,:) = reshape(mean(curIdSnrOut, 1), [1,nNexc,nMet]);
    idSnrOutStd(aa,:,:)  = reshape(std(curIdSnrOut, [], 1), [1,nNexc,nMet]);
    idSnrOutMin(aa,:,:)  = reshape(min(curIdSnrOut, [], 1), [1,nNexc,nMet]);
    idSnrOutMax(aa,:,:)  = reshape(max(curIdSnrOut, [], 1), [1,nNexc,nMet]);
    idNrmseMean(aa,:,:)  = reshape(mean(curIdNrmse, 1), [1,nNexc,nMet]);
    idNrmseStd(aa,:,:)   = reshape(std(curIdNrmse, [], 1), [1,nNexc,nMet]);
    idNrmseMin(aa,:,:)   = reshape(min(curIdNrmse, [], 1), [1,nNexc,nMet]);
    idNrmseMax(aa,:,:)   = reshape(max(curIdNrmse, [], 1), [1,nNexc,nMet]);
    idNrmseNf(aa,:,:)    = reshape(curIdNrmseNf, 1, nNexc, nMet);
    
    idPhSnrOutMean(aa,:,:) = reshape(mean(curIdPhSnrOut, 1), [1,nNexc,nMet]);
    idPhSnrOutStd(aa,:,:)  = reshape(std(curIdPhSnrOut, [], 1), [1,nNexc,nMet]);
    idPhSnrOutMin(aa,:,:)  = reshape(min(curIdPhSnrOut, [], 1), [1,nNexc,nMet]);
    idPhSnrOutMax(aa,:,:)  = reshape(max(curIdPhSnrOut, [], 1), [1,nNexc,nMet]);
    idPhNrmseMean(aa,:,:)  = reshape(mean(curIdPhNrmse, 1), [1,nNexc,nMet]);
    idPhNrmseStd(aa,:,:)   = reshape(std(curIdPhNrmse, [], 1), [1,nNexc,nMet]);
    idPhNrmseMin(aa,:,:)   = reshape(min(curIdPhNrmse, [], 1), [1,nNexc,nMet]);
    idPhNrmseMax(aa,:,:)   = reshape(max(curIdPhNrmse, [], 1), [1,nNexc,nMet]);
    idPhNrmseNf(aa,:,:)    = reshape(curIdPhNrmseNf, 1, nNexc, nMet);
    
    if useImgs
        idIms(:,:,aa,:,:,:)   = reshape(curIdIms, mtx, mtx, 1, imRep, nNexc, nMet);
        idPhIms(:,:,aa,:,:,:) = reshape(curIdPhIms, mtx, mtx, 1, imRep, nNexc, nMet);
        idImsNf(:,:,aa,:,:)   = reshape(curIdImsNf, mtx, mtx, 1, nNexc, nMet);
        idPhImsNf(:,:,aa,:,:) = reshape(curIdPhImsNf, mtx, mtx, 1, nNexc, nMet);
    end
end
spSnrOutMean = reshape(spSnrOutMean, [resSize, nMet]);
spSnrOutStd  = reshape(spSnrOutStd, [resSize, nMet]);
spSnrOutMin  = reshape(spSnrOutMin, [resSize, nMet]);
spSnrOutMax  = reshape(spSnrOutMax, [resSize, nMet]);

idSnrOutMean = reshape(idSnrOutMean, [resSize, nNexc, nMet]);
idSnrOutStd  = reshape(idSnrOutStd, [resSize, nNexc, nMet]);
idSnrOutMin  = reshape(idSnrOutMin, [resSize, nNexc, nMet]);
idSnrOutMax  = reshape(idSnrOutMax, [resSize, nNexc, nMet]);

idPhSnrOutMean = reshape(idPhSnrOutMean, [resSize, nNexc, nMet]);
idPhSnrOutStd  = reshape(idPhSnrOutStd, [resSize, nNexc, nMet]);
idPhSnrOutMin  = reshape(idPhSnrOutMin, [resSize, nNexc, nMet]);
idPhSnrOutMax  = reshape(idPhSnrOutMax, [resSize, nNexc, nMet]);

spNrmseMean = reshape(spNrmseMean, [resSize, nMet]);
spNrmseStd  = reshape(spNrmseStd, [resSize, nMet]);
spNrmseMin  = reshape(spNrmseMin, [resSize, nMet]);
spNrmseMax  = reshape(spNrmseMax, [resSize, nMet]);

idNrmseMean = reshape(idNrmseMean, [resSize, nNexc, nMet]);
idNrmseStd  = reshape(idNrmseStd, [resSize, nNexc, nMet]);
idNrmseMin  = reshape(idNrmseMin, [resSize, nNexc, nMet]);
idNrmseMax  = reshape(idNrmseMax, [resSize, nNexc, nMet]);

idPhNrmseMean = reshape(idPhNrmseMean, [resSize, nNexc, nMet]);
idPhNrmseStd  = reshape(idPhNrmseStd, [resSize, nNexc, nMet]);
idPhNrmseMin  = reshape(idPhNrmseMin, [resSize, nNexc, nMet]);
idPhNrmseMax  = reshape(idPhNrmseMax, [resSize, nNexc, nMet]);

spNrmseNf    = reshape(spNrmseNf, [resSize, nMet]);
idNrmseNf    = reshape(idNrmseNf, [resSize, nNexc, nMet]);
idPhNrmseNf  = reshape(idPhNrmseNf, [resSize, nNexc, nMet]);
spNrmseNf(:,:,2:end,:)     = [];
idNrmseNf(:,:,2:end,:,:)   = [];
idPhNrmseNf(:,:,2:end,:,:) = [];

spIms     = reshape(spIms, [mtx, mtx, resSize, imRep, nMet]);
spIms     = spIms(:,:,imFa,imT2s,imNoi,:,:);
idIms     = reshape(idIms, [mtx, mtx, resSize, imRep, nNexc, nMet]);
idIms     = idIms(:,:,imFa,imT2s,imNoi,:,:,:);
idPhIms   = reshape(idPhIms, [mtx, mtx, resSize, imRep, nNexc, nMet]);
idPhIms   = idPhIms(:,:,imFa,imT2s,imNoi,:,:,:);
spImsNf   = reshape(spImsNf, [mtx, mtx, resSize, nMet]);
spImsNf   = spImsNf(:,:,imFa,imT2s,imNoi,:);
idImsNf   = reshape(idImsNf, [mtx, mtx, resSize, nNexc, nMet]);
idImsNf   = idImsNf(:,:,imFa,imT2s,imNoi,:,:);
idPhImsNf = reshape(idPhImsNf, [mtx, mtx, resSize, nNexc, nMet]);
idPhImsNf = idPhImsNf(:,:,imFa,imT2s,imNoi,:,:);
% spImsNf(:,:,:,:,2:end,:)     = [];
% idImsNf(:,:,:,:,2:end,:,:)   = [];
% idPhImsNf(:,:,:,:,2:end,:,:) = [];

toc
endTime = datestr(now, 30);
fprintf('Finished %d iterations at %s.\n', nIter, endTime)

% return

%% Save workspace
clearvars -except trajFile nMet P fa t2s noiVar nExc dte nReps nIter fHz ...
    noiF tSp tId sp* id* im* flt resSize startTime endTime
if ispc
    compString = getenv('COMPUTERNAME');
else
    compString = getenv('HOSTNAME');
    compString = strtok(compString, '.');
end
filname = sprintf('results-%s-noFilter-%s-%s.mat', ...
    mfilename, compString, startTime);
save(filname)


return
