function [nsa, cNum] = idSimPlot(simHz, nEcho, esp, simT2s, simBw, simNum, simMag)
%IDSIMPLOT Simulates and show IDEAL performance in plots depicting
%          effective number of signal averages (NSA), phasor matrix
%          condition number, and magnitude/phase frequency response
%
%   Usage: [nsa, cNum] = idSimPlot(simHz, nEcho, [esp, simT2s, simBw, simNum, simMag])
% 
%   Inputs:
%     simHz  - Vector of chemical shifts (Hz)
%     nEcho  - Number of IDEAL echoes
%     esp    - Echo time spacing (ms)
%       If 0 or [], the echo spacing <1 ms that minimizes condition number is chosen
%       If -N, the echo spacing <1 ms that maximizes NSA for simHz(N) is chosen
%     simT2s - T2* relaxation time constant, default inf (ms)
%     simBw  - Full bandwidth to plot, default 2e3 (Hz)
%     simNum - Number of echo spacings and spectral points to simulate, default 2001
%     simMag - Vector of relative chemical shift magnitudes, default ones
%       
%   Output:
%     nsa  - Vector of effective NSA values for individual chemical shifts
%     cNum - Condition number of IDEAL phasor matrix
%
%   Literature:
%     Brodsky EK, et al. MRM 2008. DOI 10.1002/mrm.21580
%     Brodsky EK, et al. JMRI 2010. DOI 10.1002/jmri.22308
%
%   See also IDDECOMP, SHAPEBRUKERHPMR
%
% 10/2019, Keith Michel

if nargin<2, help(mfilename); return; end

%% Parse inputs
simHz     = simHz(:).';
nShifts   = numel(simHz);
if nargin<3,        esp = []; end
if isempty(esp),    esp = 0; end
if nargin<4,        simT2s = []; end
if isempty(simT2s), simT2s = inf; end
if nargin<5,        simBw = []; end
if isempty(simBw),  simBw = 2e3; end
if nargin<6,        simNum = []; end
if isempty(simNum), simNum = 2001; end
if nargin<7,        simMag = []; end
if isempty(simMag), simMag = ones(1,nShifts); end
simT2s = simT2s / 1e3; % ms -> s

%% Simulate effective NSA and condition number for range of ESP values
i2p     = 1i * 2 * pi;
simCond = zeros(1, simNum);
simNsa  = zeros(nShifts, simNum);
if esp > 0
    simEsp = [linspace(1e-3, 2*esp, simNum), esp];
    simCond(end+1)  = 0;
    simNsa(:,end+1) = 0;
else
    simEsp = linspace(1e-3, 1, simNum);
end
for ii = 1:numel(simEsp)
    simEchoTimes = (0:nEcho-1).' * simEsp(ii) / 1e3; % ms -> s
    simPhasorMatrix = exp( i2p * simEchoTimes * simHz ) .* ...
        repmat( exp( -simEchoTimes / simT2s), 1, nShifts);
    simCond(ii) = cond(simPhasorMatrix);
    simPhasorMatrixInv = pinv(simPhasorMatrix' * simPhasorMatrix);
    for jj = 1:nShifts
        simNsa(jj,ii) = 1/real(simPhasorMatrixInv(jj,jj));
    end
end
if esp > 0
    nsa       = simNsa(:,end).';
    cNum      = simCond(end);
    simEsp    = simEsp(1:end-1);
    simNsa    = simNsa(:,1:end-1);
    simCond   = simCond(1:end-1);
elseif esp == 0
    [~,idx]   = min(simCond);
    esp       = simEsp(idx);
    nsa       = simNsa(:,idx).';
    cNum      = simCond(idx);
    fprintf(1, 'Echo spacing of %g chosen to minimize condition number\n', esp);
elseif esp < 0
    [~,idx]   = max(simNsa(-esp,:));
    fprintf(1, 'Echo spacing of %g chosen to maximize NSA for %+g Hz resonance\n', ...
        simEsp(idx), simHz(-esp));
    esp  = simEsp(idx);
    nsa  = simNsa(:,idx).';
    cNum = simCond(idx);
end
echoTimes = esp * (0:nEcho-1).' / 1e3;  % ms -> s
    
figSim = figure('name','IDEAL Simulation Plots',...
    'position',[50 50 1000 900]);
movegui(figSim);

%% Effective NSA plot
subplot 321
plot(simEsp, simNsa, 'linewidth', 1.2);
hold on
line(esp*[1 1], [0 max(nsa)], 'color', 'black', ...
    'marker', 'o', 'linestyle', '--', 'linewidth', 1.2);
xlim([0 simEsp(end)])
ylim([0 nEcho+1])
xlabel('Echo Spacing (ms)', 'fontsize', 14);
set(gca, 'fontsize', 14);
grid on
title('Effective NSA', 'fontsize', 16);
    
%% Condition number plot
subplot 322
plot(simEsp, simCond, 'k');
hold on
line(esp*[1 1], [0 cNum], 'color', 'black', ...
    'marker', 'o', 'linestyle', '--', 'linewidth', 1.2);
xlim([0 simEsp(end)])
ylim([0 10])
xlabel('Echo Spacing (ms)', 'fontsize', 14);
set(gca, 'fontsize', 14);
grid on
title(sprintf('Condition Number = %.3g', cNum), ...
    'fontsize',16);
phasorMatrix = exp( i2p * echoTimes * simHz ) .* ...
    repmat( exp( echoTimes / simT2s), 1, nShifts);
simFreq = linspace(-simBw/2, simBw/2, simNum); % Hz
simDat  = exp(i2p * echoTimes * simFreq).' * pinv(phasorMatrix).';
simAbs  = abs(simDat);
simAng  = angle(simDat);

%% Frequency response magnitude plot
subplot(3,2,3:4)
hh = plot(simFreq, simAbs, 'linewidth', 1.2);
hold on
if isinf(simT2s),   simT2s = 0.02; end  % Just to make a nicer looking plot
simSpec = nscale(genspec([simMag; simHz; simT2s*ones(1, nShifts)], simBw, simNum));
plot(simFreq, simSpec, 'k', 'linewidth', 1.2);
ylabel('Magnitude')
yl = ylim;
line(repmat(simHz.', 1, 2).', repmat(yl, nShifts, 1).', ...
    'linewidth', 1.4, 'linestyle', '--');
xlim(simBw * [-1/2 1/2])
ylim(yl)
grid on
set(gca, 'fontsize', 14, 'xdir', 'reverse')
legendCell = cell(1,jj);
for jj = 1:nShifts
    legendCell{jj} = sprintf('%+.1f Hz (%.2g NSA)', simHz(jj), nsa(jj));
end
if nShifts < 7
    lh = legend(hh, legendCell, 'location', 'northoutside', 'orientation', 'horizontal', ...
        'FontSize', 9);
    ph = lh.Position + [0 0.03 0 0];
    lh.Position = [(1-ph(3))/2 ph(2:4)];
end

%% Frequency response phase plot
subplot(3,2,5:6)
plot(simFreq, simAng, '.', 'markersize', 5)
ylabel('Phase (rad)')
xlabel('Frequency (Hz)')
hold on
line(repmat(simHz.', 1, 2).', repmat(ylim, nShifts, 1).', ...
    'linewidth', 1.4, 'linestyle', '--')
xlim(simBw * [-1/2 1/2])
ylim([-3.3, 3.3])
grid on
set(gca, 'fontsize', 14, 'xdir', 'reverse', 'YTick', [-pi, 0, pi], ...
    'YTickLabel', {'-\pi', '0', '\pi'})

end

