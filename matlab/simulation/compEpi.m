% This script evaluates relative SNR and NRMSE for flyback and symmetric
% EPI readouts for a given FOV, matrix and set of gradient parameters
% (Figure 4)

clear variables
close all

%%
fov  = 4; % [cm]
mtx  = 16;
dw   = 16:16:256; % [us]
t2s  = 5:5:50;    % [ms]
opts = struct('fov', fov, ...
              'mtx', mtx, ...
              'Gmax', 600, ...  % [mT/m]
              'Smax', 4e3, ...  % [T/m/s]
              'tGrad', 8,  ...  % [us]
              'gmr', 10.71, ... % [MHz/T]
              'showPlot', false);
snrFbk = zeros(numel(dw), numel(t2s));
snrSym = zeros(numel(dw), numel(t2s));
phantm = cphantom(mtx);
k      = fftdim(phantm);
nreFbk = zeros(numel(dw), numel(t2s));
nreSym = zeros(numel(dw), numel(t2s));
for ii = 1:numel(dw)
    for jj = 1:numel(t2s)
        opts.t2 = t2s(jj);
        opts.dw = dw(ii);
        opts.epiType = false;
        [~,trajFbk(ii,jj)] = create2dEpi(opts);
        opts.epiType = true;
        [~,trajSym(ii,jj)] = create2dEpi(opts);
        snrFbk(ii,jj) = sqrt(dw(ii)) * exp(-trajFbk(ii,jj).TE ./ t2s(jj));
        snrSym(ii,jj) = sqrt(dw(ii)) * exp(-trajSym(ii,jj).TE ./ t2s(jj));
        
        kFbk = reshape(exp(-trajFbk(ii,jj).T(trajFbk(ii,jj).idxAcq) / ...
            t2s(jj)), mtx, mtx);
        pFbk = abs(ifftdim(k .* kFbk));
%         pFbk = nscale(pFbk);
        kSym = reshape(exp(-trajSym(ii,jj).T(trajSym(ii,jj).idxAcq) / ...
            t2s(jj)), mtx, mtx);
        kSym(2:2:end,:) = flipud(kSym(2:2:end,:));
        pSym = abs(ifftdim(k .* kSym));
%         pSym = nscale(pSym);
        nreFbk(ii,jj) = nrmse(pFbk, phantm);
        nreSym(ii,jj) = nrmse(pSym, phantm);
    end
end

%%
fig = figure('position', [100 50 1200 900]);
cm1 = brewermap(8, 'YlGnBu');
cm2 = brewermap(10, '*YlOrRd');
[x,y]    = meshgrid(t2s, dw);
[xx,yy]  = meshgrid(t2s(1):0.5:t2s(end), dw(1):1:dw(end));
snrFbkIm = interp2(x, y, snrFbk, xx, yy);
snrSymIm = interp2(x, y, snrSym, xx, yy);
nreFbkIm = interp2(x, y, nreFbk, xx, yy);
nreSymIm = interp2(x, y, nreSym, xx, yy);

subplot 221
contourf(yy(1:end,1), xx(1,1:end), snrFbkIm.')
ylabel('T2* (ms)')
xlabel('Dwell Time (\mus)')
colormap(gca, cm1)
colorbar
caxis([0 8])
set(gca, 'fontsize', 18, 'ydir', 'normal')
title({sprintf('Flyback EPI'), ... - %d cm FOV, %d x %d Matrix', fov, mtx, mtx), ...
    'Relative SNR'}, 'fontsize', 20)

subplot 222
contourf(yy(1:end,1), xx(1,1:end), snrSymIm.')
ylabel('T2* (ms)')
xlabel('Dwell Time (\mus)')
colormap(gca, cm1)
colorbar
caxis([0 8])
set(gca, 'fontsize', 18, 'ydir', 'normal')
title({sprintf('Symmetric EPI'), ... - %d cm FOV, %d x %d Matrix', fov, mtx, mtx), ...
    'Relative SNR'}, 'fontsize', 20)

subplot 223
contourf(yy(1:end,1), xx(1,1:end), nreFbkIm.')
ylabel('T2* (ms)')
xlabel('Dwell Time (\mus)')
colormap(gca, cm2)
colorbar
caxis([0 1])
set(gca, 'fontsize', 18, 'ydir', 'normal')
title('Normalized RMS Error', 'fontsize', 20)

subplot 224
contourf(yy(1:end,1), xx(1,1:end), nreSymIm.')
ylabel('T2* (ms)')
xlabel('Dwell Time (\mus)')
colormap(gca, cm2)
colorbar
caxis([0 1])
set(gca, 'fontsize', 18, 'ydir', 'normal')
title('Normalized RMS Error', 'fontsize', 20)
