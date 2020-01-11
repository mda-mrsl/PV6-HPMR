% Demonstration of single-band hyperpolarized C-13 pulse designs for metabolite specific MRI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectral-Spatial RF Pulse Design for MRI and MRSI MATLAB Package
%
% Authors: Adam B. Kerr and Peder E. Z. Larson
%
% (c)2007-2014 Board of Trustees, Leland Stanford Junior University and
%	The Regents of the University of California.
% All Rights Reserved.
%
% Please see the Copyright_Information and README files included with this
% package.  All works derived from this package must be properly cited.
%
%
% KAM modified 7/30/2019
%  Execute from the 'examples' directory of the Spectral-Spatial toolbox:
%   https://github.com/LarsonLab/Spectral-Spatial-RF-Pulse-Design
%  Select pulse #8 to obtain pulse shown in Figure 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add 'spectralspatial' and 'rf_tools' folders to path (mex files must be compiled)
addpath ..
addpath ../rf_tools

%% single band for high-field, small animal system
fprintf(1, '\nHere is a similar pulse design, but for a 7T small-animal system\n\n');

clear all; ss_opt([]);	ss_globals;			% Reset all options

SS_TS = 8e-6;   % KAM - 8us sampling for Bruker system

% GENERAL PULSE PARAMETERS
% ss_type = 'EP Whole';
ss_type = 'Flyback Whole';
ptype = 'ex';  % excitation pulse
opt = ss_opt({'Nucleus', 'Carbon', ...
	      'Max Duration', 20e-3, ...
	      'Max Grad', 60, ...  % G/cm       % KAM BGA12SHP High Power
	      'Max Slew', 450, ... % G/cm/ms    % KAM BGA12SHP High Power
          'Verse Fraction', 0.5, ...
          'Spect Correct', 1});

% SPECTRAL PULSE PARAMETERS  - large pass/stop bands chosen for wide
% supression regions
B0 = 7.05e4; % G
df = 0.5e-6 * B0 * SS_GAMMA; % 0.5 ppm = gamma_C13 * B0 * 0.5e-6
ea = 40; % excitation angle. If ss_design struggles iterating, try lowering
% metabolite			frequency (Hz)		freq bandwidth (Hz)		flip angle (deg)
mets(1).name = 'sb'; 	mets(1).f = 0;      mets(1).df = 2*df; 		mets(1).ang = ea;
mets(2).name = 'h1'; 	mets(2).f = 400; 	mets(2).df = 2*df; 		mets(2).ang = 0;
mets(3).name = 'l1'; 	mets(3).f = -400; 	mets(3).df = 2*df; 		mets(3).ang = 0;
mets(4).name = 'h2'; 	mets(4).f = 900; 	mets(4).df = 2*df; 		mets(4).ang = 0;
mets(5).name = 'l2'; 	mets(5).f = -900; 	mets(5).df = 2*df; 		mets(5).ang = 0;

% create vectors of angles, ripples, and band edges for input to pulse design
[fspec, a_angs, d] = create_freq_specs(mets);
fctr = 0;  % force pulse design to optimize for center of frequency specification
s_ftype = 'min';  % minimum-phase spectral filter

% SPATIAL PULSE PARAMETERS
z_thk = 1;  % thickness (cm)
z_tb = 3.5; % time-bandwidth, proportional to profile sharpness
z_ftype='ls';  % least-squares filter design
z_d1 = 0.01;  z_d2 = 0.01;  % slice profile pass and stop-band ripples, respectively

%% DESIGN THE PULSE!
% !!! Pulse #8 chosen for use on MDA 7T:
%   Fs: 1358.7 B1: 0.198G Power: 1.467e-05 G^2 ms Dur:  6.6ms
% Isolution = 8;
% [g,rf,fs,z,f,mxy] = ...
%     ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
% 	      z_ftype, s_ftype, ss_type, fctr, [], [], Isolution);
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr, [], []);
close(gcf)
% Rescale B1 to 90 deg excitation, RF and gradient waveforms to preferred units
rf       = rf * 90/ea;
grad_mTm = g * 10;                              % G/cm -> mT/m
slew_Tms = abs(diff(grad_mTm)) / SS_TS / 1e3;   % T/m/s
rf_uT    = rf * 100;                            % G -> uT

%% KAM custom plots for HP pyr imaging
tax    = (1:numel(rf)) * SS_TS * 1e3;
rf_rad = angle(rf);
rfn    = abs(rf) / max(abs(rf(:)));
met_name = {'bic', 'urea', 'pyr', 'ala', 'pyh', 'lac'};
met_ppm  = [161, 164, 171, 177, 179.5, 183.2];
met_Hz   = (met_ppm - 171) * SS_GAMMA * B0 * 1e-6;
met_mag  = [2, 2, 10, 2, 3, 6];
colors   = lines(numel(met_Hz));
bw = 4e3;
[fp,zp,mp] = ss_plot(g, rf, SS_TS, ptype, z_thk*2, [-bw/2, bw/2], ...
    SS_GAMMA, met_Hz); % redo plot
fg1 = gcf;
fg1.Position = [200 50 1200 950];
fg2 = figure('position', [100 100 1500 300]);
fg3 = figure('position', [50 50 600 850]);
rfp = sum(abs(mp),1);
rfp = rfp - min(rfp);
rfp = rfp / max(rfp);
rfpm = max(abs(mp));
rfpm = rfpm - min(rfpm);
rfpm = rfpm / max(rfpm);
spec = abs(genspec([met_mag; met_Hz; 0.01*ones(1,6)], 4e3, 1e3));
spec = spec - min(spec);
spec = spec / max(spec);
plot(fp, rfp, 'k', 'linewidth', 1.5)
hold on
plot(fp, rfpm, '--k', 'linewidth', 1.5)
excFrac = zeros(numel(met_Hz)*[1 1]);
fpp = linspace(-bw/2, bw/2, 1e3);
rfp = interp1(fp, rfp, fpp);
spp = zeros(1, 1e3);
for ii = 1:numel(met_Hz)
    spp(find(fpp>met_Hz(ii), 1)) = met_mag(ii)/10;
end
for ii = 1:numel(met_Hz)
    excFrac(ii,:) = interp1(fpp, rfp, met_Hz - met_Hz(ii));

    lift = (ii-1)/numel(met_Hz);
    offs = numel(fpp)/2 - find(fpp>met_Hz(ii), 1);
    curs = circshift(spec, offs);
    curd = circshift(spp, offs);

    axes(fg1.Children(ii))
    cla
    area(fpp, curs, 'linewidth', 1, 'facealpha', 1, ...
        'edgecolor', 'none', 'facecolor', colors(ii,:));
    hold on
    plot(fpp, rfp, 'k');
    xlabel('Frequency (Hz)')
    set(gca, 'xdir', 'reverse', 'ytick', [])
    title(met_name{ii})
    figure(fg2)
    subplot(1, numel(met_Hz), ii)
    plot(fpp, rfp, 'k', 'linewidth', 1.5);
    hold on
    plot(fpp, curs, 'linewidth', 1.5);
    xlabel('Frequency (Hz)')
    set(gca, 'fontsize', 15, 'xdir', 'reverse', 'ytick', [])
    title(met_name{ii}, 'fontsize', 16)
    figure(fg3)
    plot(fpp, lift+curs/numel(met_Hz), 'linewidth', 1.5);
end
figure(fg3)
xlabel('Frequency (Hz)')
grid on
set(gca, 'fontsize', 15, 'xdir', 'reverse', 'ytick', [])
axes(fg1.Children(numel(met_Hz)+4))
cla
plot(tax, grad_mTm, 'linewidth', 1);
xlabel('Time (ms)')
ylabel('G (mT/m)')
xlim([0 ceil(round(tax(end)))])
grid on
title('Gradient')
fg1.Children(1).Position = [0.1300 0.58 0.7750 0.14];
axes(fg1.Children(numel(met_Hz)+5))
cla
plot(tax, real(rf_uT), 'linewidth', 1);
hold on
plot(tax, imag(rf_uT), 'linewidth', 1);
plot(tax, abs(rf_uT), '--', 'linewidth', 1);
xlabel('Time (ms)')
ylabel('B1 (uT)')
xlim([0 ceil(round(tax(end)))])
grid on
title('RF Pulse')
fg1.Children(1).Position = [0.1300 0.8 0.7750 0.14];

% return

%%
fg4 = figure('position', [200 50 1200 900]);
met_name = {'lactate', {'pyruvate', 'hydrate'}, 'alanine', 'pyruvate', 'urea', 'bicarb'};
met_ppm  = [183.2, 179.5, 177, 171, 164, 161];
met_Hz   = (met_ppm - 171) * SS_GAMMA * B0 * 1e-6;
met_mag  = [6, 3, 2, 10, 2, 2];
colors   = lines(6);
spec = abs(genspec([met_mag; met_Hz; 0.01*ones(1,6)], 4e3, 1e3));
spec = spec - min(spec);
spec = spec / max(spec);
subplot(2,3,2:3)
plot(tax, real(rf_uT), 'linewidth', 1);
hold on
plot(tax, imag(rf_uT), 'linewidth', 1);
plot(tax, abs(rf_uT), '--', 'linewidth', 1);
ylabel('B1 (uT)')
xlim([0 ceil(round(tax(end)))])
grid on
legend('Real', 'Imaginary', 'Magnitude')
set(gca, 'fontsize', 14)
title('RF Pulse', 'fontsize', 15)
subplot(2,3,5:6)
plot(tax, grad_mTm, 'linewidth', 1);
xlabel('Time (ms)')
ylabel('G (mT/m)')
xlim([0 ceil(round(tax(end)))])
grid on
set(gca, 'fontsize', 14)
title('Gradient', 'fontsize', 15)
subplot(2,3,[1,4])
plot(fpp, rfp, 'k', 'linewidth', 1.5)
hold on
xx = [1075, 1075, 1075, -300, -350, -320];
yy = [ 0.9, 0.80, 0.62, 0.4, 0.25, 0.1];
for ii = 1:numel(met_Hz)
    lift = (numel(met_Hz)-ii)/numel(met_Hz);
    offs = numel(fpp)/2 - find(fpp>met_Hz(ii), 1);
    curs = circshift(spec, offs);
    plot(fpp, lift+curs/numel(met_Hz), 'linewidth', 1.5);
    text(xx(ii), yy(ii), met_name{ii}, 'fontsize', 14, ...
        'color', colors(ii,:))
end
xlabel('Frequency (Hz)')
grid on
set(gca, 'fontsize', 14, 'YTick', [], 'xdir', 'reverse')


%% Write rf and gradient waveform to files
max_rfuT = 100*max(abs(rf));
bpB1uT   = 1 / (SS_GAMMA * 4e-5);       % B1 for 1 ms, 90 deg bp [uT]
powerFactor = (max_rfuT / bpB1uT)^2;    % RF power is propto B1^2
rfn = 100 * rf / max_rfuT;
sInt = sqrt( sum(real(rfn))^2 + sum(imag(rfn))^2) / numel(rfn);
pInt = rfn * rfn' / numel(rfn);
length_ms = 1e3 * SS_TS * numel(rfn);

rfn = rfn * 100;
mag = abs(rfn);
ph_deg = 180 + angle(rfn) * 180/pi;
grad_mTm = g * 10;

slOffs  = 0; % mm
nSlices = numel(slOffs);

curTime = clock;
curTime(6) = round(curTime(6));
filName = sprintf('SpSp_13C_sb_%dmm_%d%02d%02d',...
    round(10*z_thk),curTime(1:3));
waveDir = fullfile(pwd,'wave');
gpDir = fullfile(pwd,'gp');
jcampHeader = {
    sprintf('##TITLE= %s',fullfile(waveDir,filName));
    '##JCAMP-DX= 5.00 Bruker JCAMP library';
    '##DATA TYPE= Shape Data';
    '##ORIGIN= MDACC';
    '##OWNER= <keith>';
    sprintf('##DATE= %d/%02d/%02d',curTime(1:3));
    sprintf('##TIME= %02d:%02d:%02d',curTime(4:6));
    sprintf('##MINX= %.6e',min(mag));
    sprintf('##MAXX= %.6e',max(mag));
    sprintf('##MINY= %.6e',min(ph_deg));
    sprintf('##MAXY= %.6e',max(ph_deg));
    '##$SHAPE_EXMODE= Excitation';
    '##$SHAPE_TOTROT= 9.000000e+01';
    '##$SHAPE_BWFAC= ';
    sprintf('##$SHAPE_INTEGFAC= %.8e',sInt);
    '##$SHAPE_REPHFAC=50';
    '##$SHAPE_TYPE=conventional';
    '##$SHAPE_MODE= 0';
    sprintf('##MAXB1 = %.5e uT',max_rfuT);
    sprintf('##NPOINTS= %d',numel(rf));
    sprintf('##DURATION= %.5e ms',length_ms);
    sprintf('##NUCLEUS= %s',opt{1,2});
    sprintf('##FIELD= %.5e T',B0*1e-4);
    sprintf('##MAXGRAD= %.5e mT/m',opt{2,2}*10);
    sprintf('##MAXSLEW= %.5e T/m/s',opt{3,2}*10);
    sprintf('##SLICEWIDTH= %.5e cm',z_thk);
    sprintf('##SLICEOFFSET= %.5e mm',0);
    '##XYPOINTS= (XY..XY)'};

if ~exist(waveDir,'dir'), mkdir(waveDir); end
rfFile = fopen(fullfile(waveDir,sprintf('%s.exc',filName)),'w');
for ii=1:numel(jcampHeader)
    fprintf(rfFile,'%s\n',jcampHeader{ii});
end
for ii=1:numel(rf)
    fprintf(rfFile,'%.6e, %.6e\n',mag(ii),ph_deg(ii));
end
fprintf(rfFile,'##END= \n');
fclose(rfFile);

metaFile = fopen(fullfile(waveDir,sprintf('%s.meta',filName)),'w');
fprintf(metaFile,'%.6e %.6e %d', ...
    length_ms, sInt, numel(rf));
fclose(metaFile);

if ~exist(gpDir,'dir'), mkdir(gpDir); end
gradFile = fopen(fullfile(gpDir,sprintf('%s.gp',filName)),'w');
for ii=1:numel(g)
    fprintf(gradFile,'%.6e\n',grad_mTm(ii));
end
fclose(gradFile);

save(fullfile(waveDir,sprintf('%s.mat',filName)));
return

