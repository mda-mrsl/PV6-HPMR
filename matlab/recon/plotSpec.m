function plotSpec( spec, bw, logScale )
%PLOTSPEC plotting for spectroscopic MR data
%
%   Usage:  plotSpec(spec, bw, logScale)
%
%       where spec is an array containing spectroscopic data
%               if complex, magnitude is used
%               if >1 dimension, all columns are shown on same axis
%             bw is the full bandwidth in Hz
%               if omitted, x axis is unlabeled
%             logScale is a boolean indicating whether to use semilogy
%               if omitted, default is false
%
%   See also PLOT, SEMILOGY, SHOWMONTAGE
%
%   06/2019, Keith Michel

%%
if ~nargin,     help(mfilename); return, end
if nargin<2,    bw = []; end
if nargin<3,    logScale = false; end

%%
n    = size(spec, 1);
spec = reshape(spec, n, []);
if ~isreal(spec)
    spec = abs(spec);
end

if isempty(bw)
    fax = fliplr(1:n);
else
    fax = linspace(bw/2, -bw/2, n);
end

figure
set(gca, 'colororder', winter(size(spec, 2)), 'nextplot', 'replacechildren')
if logScale
    semilogy(fax, spec);
else
    plot(fax, spec);
end
if ~isempty(bw)
    xlabel('Frequency (Hz)')
end
set(gca, 'xdir', 'reverse')
grid on
shg
end
