% This script plots relative SNR and NRMSE results for spectral encoding 
% via spectral-spatial excitation and IDEAL encoding (Figure 5).
% 
% Some plots require the boundedline function:
% https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m

clear variables
close all

% The .mat file containing results to plot:
matFile = 'results-compIdSpsp-noFilter-pigpen-20190919T184558';
load(matFile)

%% Contour plots of SNR and NRMSE averaged across metabs
figure('position', [100 100 1200 900])
cm1 = brewermap(64, 'YlGnBu');
cm2 = brewermap(64, '*YlOrRd');
[~,rNvar]  = min(abs(noiVar - 0.05));
[~,rNexc]  = min(abs(nExc - 8));
[x,y]      = meshgrid(t2s, fa);
[xx,yy]    = meshgrid(t2s(1):0.5:t2s(end), fa(1):0.5:fa(end));
snrSpCont  = interp2(x, y, mean(spSnrOutMean(:,:,rNvar,:), 4), xx, yy);
snrIdCont  = interp2(x, y, mean(idSnrOutMean(:,:,rNvar,rNexc,:), 5), xx, yy);
nreSpCont  = interp2(x, y, mean(spNrmseMean(:,:,rNvar,:), 4), xx, yy);
nreIdCont  = interp2(x, y, mean(idNrmseMean(:,:,rNvar,rNexc,:), 5), xx, yy);

subplot 221
contourf(yy(1:end,1), xx(1,1:end), snrSpCont.')
xlabel('Effective Excitation Angle (deg)')
ylabel('T2* (ms)')
colormap(gca, cm1)
colorbar
set(gca, 'fontsize', 18, 'ydir', 'normal')
title({'Spectral Spatial Excitation', 'SNR'}, 'fontsize', 20)

subplot 222
contourf(yy(1:end,1), xx(1,1:end), snrIdCont.')
xlabel('Effective Excitation Angle (deg)')
ylabel('T2* (ms)')
colormap(gca, cm1)
colorbar
set(gca, 'fontsize', 18, 'ydir', 'normal')
title({'Multi Echo IDEAL', 'SNR'}, 'fontsize', 20)

subplot 223
contourf(yy(1:end,1), xx(1,1:end), nreSpCont.')
xlabel('Effective Excitation Angle (deg)')
ylabel('T2* (ms)')
colormap(gca, cm2)
colorbar
set(gca, 'fontsize', 18, 'ydir', 'normal')
title('Normalized RMS Error', 'fontsize', 20)

subplot 224
contourf(yy(1:end,1), xx(1,1:end), nreIdCont.')
xlabel('Effective Excitation Angle (deg)')
ylabel('T2* (ms)')
colormap(gca, cm2)
colorbar
set(gca, 'fontsize', 18, 'ydir', 'normal')
title('Normalized RMS Error', 'fontsize', 20)

%% Plot SNR and NRMSE results - boundedline plots averaged across metabs
% Standard values: fa=20deg, t2s=20ms, noiVar=0.05, nExc=8
addpath('boundedline')
[~,rFa]    = min(abs(fa - 20));
[~,rT2s]   = min(abs(t2s - 20));
[~,rNvar]  = min(abs(noiVar - 0.05));
[~,rNexc]  = min(abs(nExc - 8));
cmaps = lines(3);

% figure('position', [100 100 1600 800])
figure('position', [100 100 1200 900])

% Flip angle
% subplot 231
subplot 221
meanSp = mean(squeeze(spSnrOutMean(:,rT2s,rNvar,:)), 2);
stdSp  = mean(squeeze(spSnrOutStd(:,rT2s,rNvar,:)), 2);
% minSp  = min(squeeze(spSnrOutMin(:,rT2s,rNvar,:)), [], 2);
% maxSp  = max(squeeze(spSnrOutMax(:,rT2s,rNvar,:)), [], 2);
meanId = mean(squeeze(idSnrOutMean(:,rT2s,rNvar,rNexc,:)), 2);
stdId  = mean(squeeze(idSnrOutStd(:,rT2s,rNvar,rNexc,:)), 2);
% minId  = min(squeeze(idSnrOutMin(:,rT2s,rNvar,rNexc,:)), [], 2);
% maxId  = max(squeeze(idSnrOutMax(:,rT2s,rNvar,rNexc,:)), [], 2);
meanIdPh = mean(squeeze(idPhSnrOutMean(:,rT2s,rNvar,rNexc,:)), 2);
stdIdPh  = mean(squeeze(idPhSnrOutStd(:,rT2s,rNvar,rNexc,:)), 2);
% minIdPh  = min(squeeze(idPhSnrOutMin(:,rT2s,rNvar,rNexc,:)), [], 2);
% maxIdPh  = max(squeeze(idPhSnrOutMax(:,rT2s,rNvar,rNexc,:)), [], 2);
center = [meanSp, meanId, meanIdPh].';
width = cat(3, stdSp, stdId, stdIdPh);
center = [meanSp, meanId].';
width = cat(3, stdSp, stdId);
h = boundedline(fa, center, width, 'alpha', 'cmap', cmaps);
for jj = 1:numel(h)
    h(jj).LineWidth = 1.5;
end
h(2).LineStyle = '--';
% l = legend(h, 'Spectral Excitation', 'Multi-Echo IDEAL', 'Multi-Echo IDEAL Thermal', ...
l = legend(h, 'Spectral Excitation', 'Multi-Echo IDEAL', ...
    'location', 'best', 'autoupdate', 'off');
xlabel('Effective Excitation Angle (deg)')
ylabel('SNR')
xlim([0, 90])
% ylim([0, 25])
line(fa(rFa)*[1 1], ylim, 'linewidth', 1.5, 'linestyle', '--', 'color', 0.5*ones(1,3))
set(gca, 'fontsize', 15, 'XTick', 0:30:90)
l.FontSize = 14;
set(l.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1;1;1;0.5]))
grid on

% T2*
% subplot 232
subplot 222
meanSp = mean(squeeze(spSnrOutMean(rFa,:,rNvar,:)), 2);
stdSp  = mean(squeeze(spSnrOutStd(rFa,:,rNvar,:)), 2);
% minSp  = min(squeeze(spSnrOutMin(rFa,:,rNvar,:)), [], 2);
% maxSp  = max(squeeze(spSnrOutMax(rFa,:,rNvar,:)), [], 2);
meanId = mean(squeeze(idSnrOutMean(rFa,:,rNvar,rNexc,:)), 2);
stdId  = mean(squeeze(idSnrOutStd(rFa,:,rNvar,rNexc,:)), 2);
% minId  = min(squeeze(idSnrOutMin(rFa,:,rNvar,rNexc,:)), [], 2);
% maxId  = max(squeeze(idSnrOutMax(rFa,:,rNvar,rNexc,:)), [], 2);
meanIdPh = mean(squeeze(idPhSnrOutMean(rFa,:,rNvar,rNexc,:)), 2);
stdIdPh  = mean(squeeze(idPhSnrOutStd(rFa,:,rNvar,rNexc,:)), 2);
% minIdPh  = min(squeeze(idPhSnrOutMin(rFa,:,rNvar,rNexc,:)), [], 2);
% maxIdPh  = max(squeeze(idPhSnrOutMax(rFa,:,rNvar,rNexc,:)), [], 2);
% center = [meanSp, meanId, meanIdPh].';
% width  = abs(cat(3, [minSp, maxSp] - meanSp, ...
%     [minId, maxId] - meanId, ...
%     [minIdPh, maxIdPh] - meanIdPh));
% h = boundedline(t2s, center, width, 'alpha');
% for jj = 1:numel(h)
%     h(jj).LineWidth = 1.5;
% end
% hold on
% width = cat(3, stdSp, stdId, stdIdPh);
% boundedline(t2s, center, width, 'alpha', ...
%     'transparency', 0.5, 'cmap', cmaps);
% legend(h, 'Spectral-Spatial', 'IDEAL - HP', 'IDEAL - Thermal', ...
%     'location', 'south', 'autoupdate', 'off')
center = [meanSp, meanId, meanIdPh].';
width = cat(3, stdSp, stdId, stdIdPh);
center = [meanSp, meanId].';
width = cat(3, stdSp, stdId);
h = boundedline(t2s, center, width, 'alpha', 'cmap', cmaps);
for jj = 1:numel(h)
    h(jj).LineWidth = 1.5;
end
h(2).LineStyle = '--';
% l = legend(h, 'Spectral Excitation', 'Multi-Echo IDEAL', 'Multi-Echo IDEAL Thermal', ...
l = legend(h, 'Spectral Excitation', 'Multi-Echo IDEAL', ...
    'location', 'best', 'autoupdate', 'off');
xlabel('T2* Relaxation Time (ms)')
ylabel('SNR')
xlim([0, 100])
% ylim([0, 15])
line(t2s(rT2s)*[1 1], ylim, 'linewidth', 1.5, 'linestyle', '--', 'color', 0.5*ones(1,3))
set(gca, 'fontsize', 15)
l.FontSize = 14;
set(l.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1;1;1;0.5]))
grid on

% Noise Variance
% subplot 233
subplot 223
meanSp = mean(squeeze(spSnrOutMean(rFa,rT2s,:,:)), 2);
stdSp  = mean(squeeze(spSnrOutStd(rFa,rT2s,:,:)), 2);
% minSp  = min(squeeze(spSnrOutMin(rFa,rT2s,:,:)), [], 2);
% maxSp  = max(squeeze(spSnrOutMax(rFa,rT2s,:,:)), [], 2);
meanId = mean(squeeze(idSnrOutMean(rFa,rT2s,:,rNexc,:)), 2);
stdId  = mean(squeeze(idSnrOutStd(rFa,rT2s,:,rNexc,:)), 2);
% minId  = min(squeeze(idSnrOutMin(rFa,rT2s,:,rNexc,:)), [], 2);
% maxId  = max(squeeze(idSnrOutMax(rFa,rT2s,:,rNexc,:)), [], 2);
meanIdPh = mean(squeeze(idPhSnrOutMean(rFa,rT2s,:,rNexc,:)), 2);
stdIdPh  = mean(squeeze(idPhSnrOutStd(rFa,rT2s,:,rNexc,:)), 2);
% minIdPh  = min(squeeze(idPhSnrOutMin(rFa,rT2s,:,rNexc,:)), [], 2);
% maxIdPh  = max(squeeze(idPhSnrOutMax(rFa,rT2s,:,rNexc,:)), [], 2);
% center = [meanSp, meanId, meanIdPh].';
% width  = abs(cat(3, [minSp, maxSp] - meanSp, ...
%     [minId, maxId] - meanId, ...
%     [minIdPh, maxIdPh] - meanIdPh));
% h = boundedline(noiVar, center, width, 'alpha');
% for jj = 1:numel(h)
%     h(jj).LineWidth = 1.5;
% end
% hold on
% width = cat(3, stdSp, stdId, stdIdPh);
% boundedline(noiVar, center, width, 'alpha', ...
%     'transparency', 0.5, 'cmap', cmaps);
% legend(h, 'Spectral-Spatial', 'IDEAL - HP', 'IDEAL - Thermal', ...
%     'location', 'northwest', 'autoupdate', 'off')
center = [meanSp, meanId, meanIdPh].';
width = cat(3, stdSp, stdId, stdIdPh);
center = [meanSp, meanId].';
width = cat(3, stdSp, stdId);
h = boundedline(noiVar, center, width, 'alpha', 'cmap', cmaps);
for jj = 1:numel(h)
    h(jj).LineWidth = 1.5;
end
h(2).LineStyle = '--';
% l = legend(h, 'Spectral Excitation', 'Multi-Echo IDEAL', 'Multi-Echo IDEAL Thermal', ...
l = legend(h, 'Spectral Excitation', 'Multi-Echo IDEAL', ...
    'location', 'best', 'autoupdate', 'off');
xlabel('Noise Variance')
ylabel('SNR')
xlim([0, 0.5])
% ylim([0, 20])
line(noiVar(rNvar)*[1 1], ylim, 'linewidth', 1.5, 'linestyle', '--', 'color', 0.5*ones(1,3))
set(gca, 'fontsize', 15)
l.FontSize = 14;
set(l.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1;1;1;0.5]))
grid on

% T2* noise-free
% subplot 235
subplot 224
meanSp = mean(squeeze(spNrmseNf(rFa,:,1,:)), 2);
meanId = mean(squeeze(idNrmseNf(rFa,:,1,rNexc,:)), 2);
meanIdPh = mean(squeeze(idPhNrmseNf(rFa,:,1,rNexc,:)), 2);
% center = [meanSp, meanId, meanIdPh].';
% h = plot(t2s, center, 'linewidth', 1.5);
% legend(h, 'Spectral-Spatial', 'IDEAL - HP', 'IDEAL - Thermal', ...
%     'location', 'south', 'autoupdate', 'off')
center = [meanSp, meanId, meanIdPh].';
center = [meanSp, meanId].';
h = plot(t2s, center, 'linewidth', 2.2);
h(2).LineStyle = '--';
% h(3).LineStyle = '-.';
% l = legend(h, 'Spectral Excitation', 'Multi-Echo IDEAL', 'Multi-Echo IDEAL Thermal', ...
l = legend(h, 'Spectral Excitation', 'Multi-Echo IDEAL', ...
    'location', 'best', 'autoupdate', 'off');
xlabel('T2* Relaxation Time (ms)')
ylabel({'NRMSE', '(Noise Free)'})
% ylabel('NRMSE')
xlim([0, 90])
ylim([0, 0.8])
line(t2s(rT2s)*[1 1], ylim, 'linewidth', 1.5, 'linestyle', '--', 'color', 0.5*ones(1,3))
set(gca, 'fontsize', 15)
l.FontSize = 14;
set(l.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1;1;1;0.5]))
grid on

return

% Phantom Image
subplot 236
imagesc(P)
cmap = [zeros(1,3); linspecer(6)];
colormap(gca, cmap)
axis square
title('Simulation Phantom', 'fontsize', 15)
set(gca, 'fontsize', 15)

return

%% Plot SNR results - boundedline plots for each metab
addpath('boundedline')

% Standard values: fa=20deg, t2s=20ms, noiVar=0.05, nExc=8
[~,rFa]    = min(abs(fa - 20));
[~,rT2s]   = min(abs(t2s - 20));
[~,rNvar]  = min(abs(noiVar - 0.05));
[~,rNexc]  = min(abs(nExc - 8));
cmaps  = lines(3);
mets   = {'Lactate', 'Pyruvate-Hydrate', 'Alanine', ...
    'Pyruvate', 'Urea', 'Bicarbonate'};

% Flip angle
meanSp = squeeze(spSnrOutMean(:,rT2s,rNvar,:));
stdSp  = squeeze(spSnrOutStd(:,rT2s,rNvar,:));
minSp  = squeeze(spSnrOutMin(:,rT2s,rNvar,:));
maxSp  = squeeze(spSnrOutMax(:,rT2s,rNvar,:));
meanId = squeeze(idSnrOutMean(:,rT2s,rNvar,rNexc,:));
stdId  = squeeze(idSnrOutStd(:,rT2s,rNvar,rNexc,:));
minId  = squeeze(idSnrOutMin(:,rT2s,rNvar,rNexc,:));
maxId  = squeeze(idSnrOutMax(:,rT2s,rNvar,rNexc,:));
meanIdPh = squeeze(idPhSnrOutMean(:,rT2s,rNvar,rNexc,:));
stdIdPh  = squeeze(idPhSnrOutStd(:,rT2s,rNvar,rNexc,:));
minIdPh  = squeeze(idPhSnrOutMin(:,rT2s,rNvar,rNexc,:));
maxIdPh  = squeeze(idPhSnrOutMax(:,rT2s,rNvar,rNexc,:));
figure('position', [100 100 1200 900]);
for ii = 1:nMet
    subplot(2,3,ii)
    center = [meanSp(:,ii), meanId(:,ii), meanIdPh(:,ii)].';
    width  = abs(cat(3, [minSp(:,ii), maxSp(:,ii)] - meanSp(:,ii), ...
        [minId(:,ii), maxId(:,ii)] - meanId(:,ii), ...
        [minIdPh(:,ii), maxIdPh(:,ii)] - meanIdPh(:,ii)));
    hl = boundedline(fa, center, width, 'alpha');
    for jj = 1:numel(hl)
        hl(jj).LineWidth = 1.3;
    end
    hold on
    width = cat(3, stdSp(:,ii), stdId(:,ii), stdIdPh(:,ii));
    boundedline(fa, center, width, 'alpha', ...
        'transparency', 0.5, 'cmap', cmaps);
    legend(hl, 'Spectral-Spatial', 'IDEAL - HP', 'IDEAL - Thermal', ...
        'location', 'south')
    xlabel('Effective Excitation Angle (deg)')
    ylabel('SNR')
    xlim([0, 90])
    set(gca, 'fontsize', 13, 'XTick', 0:30:90)
    grid on
    title(sprintf('%s: %d Hz', mets{ii}, fHz(ii)), 'fontsize', 15)
end
% print(gcf,'compIdSpsp-fa','-dpng');

% T2*
meanSp = squeeze(spSnrOutMean(rFa,:,rNvar,:));
stdSp  = squeeze(spSnrOutStd(rFa,:,rNvar,:));
minSp  = squeeze(spSnrOutMin(rFa,:,rNvar,:));
maxSp  = squeeze(spSnrOutMax(rFa,:,rNvar,:));
meanId = squeeze(idSnrOutMean(rFa,:,rNvar,rNexc,:));
stdId  = squeeze(idSnrOutStd(rFa,:,rNvar,rNexc,:));
minId  = squeeze(idSnrOutMin(rFa,:,rNvar,rNexc,:));
maxId  = squeeze(idSnrOutMax(rFa,:,rNvar,rNexc,:));
meanIdPh = squeeze(idPhSnrOutMean(rFa,:,rNvar,rNexc,:));
stdIdPh  = squeeze(idPhSnrOutStd(rFa,:,rNvar,rNexc,:));
minIdPh  = squeeze(idPhSnrOutMin(rFa,:,rNvar,rNexc,:));
maxIdPh  = squeeze(idPhSnrOutMax(rFa,:,rNvar,rNexc,:));
figure('position', [100 100 1200 900]);
for ii = 1:nMet
    subplot(2,3,ii)
    center = [meanSp(:,ii), meanId(:,ii), meanIdPh(:,ii)].';
    width  = abs(cat(3, [minSp(:,ii), maxSp(:,ii)] - meanSp(:,ii), ...
        [minId(:,ii), maxId(:,ii)] - meanId(:,ii), ...
        [minIdPh(:,ii), maxIdPh(:,ii)] - meanIdPh(:,ii)));
    hl = boundedline(t2s, center, width, 'alpha');
    for jj = 1:numel(hl)
        hl(jj).LineWidth = 1.3;
    end
    hold on
    width = cat(3, stdSp(:,ii), stdId(:,ii), stdIdPh(:,ii));
    boundedline(t2s, center, width, 'alpha', ...
        'transparency', 0.5, 'cmap', cmaps);
    legend(hl, 'Spectral-Spatial', 'IDEAL - HP', 'IDEAL - Thermal', ...
        'location', 'south')
    xlabel('T2* (ms)')
    ylabel('SNR')
    xlim([0, 50])
    set(gca, 'fontsize', 13)
    grid on
    title(sprintf('%s: %d Hz', mets{ii}, fHz(ii)), 'fontsize', 15)
end
% print(gcf,'compIdSpsp-t2s','-dpng');

% Noise Variance
meanSp = squeeze(spSnrOutMean(rFa,rT2s,:,:));
stdSp  = squeeze(spSnrOutStd(rFa,rT2s,:,:));
minSp  = squeeze(spSnrOutMin(rFa,rT2s,:,:));
maxSp  = squeeze(spSnrOutMax(rFa,rT2s,:,:));
meanId = squeeze(idSnrOutMean(rFa,rT2s,:,rNexc,:));
stdId  = squeeze(idSnrOutStd(rFa,rT2s,:,rNexc,:));
minId  = squeeze(idSnrOutMin(rFa,rT2s,:,rNexc,:));
maxId  = squeeze(idSnrOutMax(rFa,rT2s,:,rNexc,:));
meanIdPh = squeeze(idPhSnrOutMean(rFa,rT2s,:,rNexc,:));
stdIdPh  = squeeze(idPhSnrOutStd(rFa,rT2s,:,rNexc,:));
minIdPh  = squeeze(idPhSnrOutMin(rFa,rT2s,:,rNexc,:));
maxIdPh  = squeeze(idPhSnrOutMax(rFa,rT2s,:,rNexc,:));
figure('position', [100 100 1200 900]);
for ii = 1:nMet
    subplot(2,3,ii)
    center = [meanSp(:,ii), meanId(:,ii), meanIdPh(:,ii)].';
    width  = abs(cat(3, [minSp(:,ii), maxSp(:,ii)] - meanSp(:,ii), ...
        [minId(:,ii), maxId(:,ii)] - meanId(:,ii), ...
        [minIdPh(:,ii), maxIdPh(:,ii)] - meanIdPh(:,ii)));
    hl = boundedline(noiVar, center, width, 'alpha');
    for jj = 1:numel(hl)
        hl(jj).LineWidth = 1.3;
    end
    hold on
    width = cat(3, stdSp(:,ii), stdId(:,ii), stdIdPh(:,ii));
    boundedline(noiVar, center, width, 'alpha', ...
        'transparency', 0.5, 'cmap', cmaps);
    legend(hl, 'Spectral-Spatial', 'IDEAL - HP', 'IDEAL - Thermal', ...
        'location', 'south')
    xlabel('Noise Variance')
    ylabel('SNR')
    xlim([0, 0.5])
    set(gca, 'fontsize', 13)
    grid on
    title(sprintf('%s: %d Hz', mets{ii}, fHz(ii)), 'fontsize', 15)
end
% print(gcf,'compIdSpsp-noiVar','-dpng');

return

% Number of IDEAL excitations
meanSp = repmat(squeeze(spSnrOutMean(rFa,rT2s,rNvar,:)).', numel(nExc), 1);
stdSp  = repmat(squeeze(spSnrOutStd(rFa,rT2s,rNvar,:)).', numel(nExc), 1);
minSp  = repmat(squeeze(spSnrOutMin(rFa,rT2s,rNvar,:)).', numel(nExc), 1);
maxSp  = repmat(squeeze(spSnrOutMax(rFa,rT2s,rNvar,:)).', numel(nExc), 1);
meanId = squeeze(idSnrOutMean(rFa,rT2s,rNvar,:,:));
stdId  = squeeze(idSnrOutStd(rFa,rT2s,rNvar,:,:));
minId  = squeeze(idSnrOutMin(rFa,rT2s,rNvar,:,:));
maxId  = squeeze(idSnrOutMax(rFa,rT2s,rNvar,:,:));
meanIdPh = squeeze(idPhSnrOutMean(rFa,rT2s,rNvar,:,:));
stdIdPh  = squeeze(idPhSnrOutStd(rFa,rT2s,rNvar,:,:));
minIdPh  = squeeze(idPhSnrOutMin(rFa,rT2s,rNvar,:,:));
maxIdPh  = squeeze(idPhSnrOutMax(rFa,rT2s,rNvar,:,:));
figure('position', [100 100 1200 900]);
for ii = 1:nMet
    subplot(2,3,ii)
    center = [meanSp(:,ii), meanId(:,ii), meanIdPh(:,ii)].';
    width  = abs(cat(3, [minSp(:,ii), maxSp(:,ii)] - meanSp(:,ii), ...
        [minId(:,ii), maxId(:,ii)] - meanId(:,ii), ...
        [minIdPh(:,ii), maxIdPh(:,ii)] - meanIdPh(:,ii)));
    hl = boundedline(nExc, center, width, 'alpha');
    for jj = 1:numel(hl)
        hl(jj).LineWidth = 1.3;
    end
    hold on
    width = cat(3, stdSp(:,ii), stdId(:,ii), stdIdPh(:,ii));
    boundedline(nExc, center, width, 'alpha', ...
        'transparency', 0.5, 'cmap', cmaps);
    legend(hl, 'Spectral-Spatial', 'IDEAL - HP', 'IDEAL - Thermal', ...
        'location', 'south')
    xlabel('Number of IDEAL Excitations')
    ylabel('SNR')
    set(gca, 'fontsize', 13)
    grid on
    title(sprintf('%s: %d Hz', mets{ii}, fHz(ii)), 'fontsize', 15)
end
% print(gcf,'compIdSpsp-nExc','-dpng');
