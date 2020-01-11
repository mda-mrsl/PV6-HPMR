function s = genspec(p,bw,n)
%GENSPEC creates synthetic NMR spectrum
%
%   Usage: s = genspec(p, f, n)
%
%       where p is a 3xN matrix with parameters of N Lorentzian peaks
%               1st row specifies amplitudes
%               2nd row specifies center frequencies [same units as bw]
%               3rd row specifies T2* relaxation rate [inverse units to bw]
%             bw is a the full sampling bandwidth [same units as p(2,:)]
%             n is the number of points
%
%   07/2019, Keith Michel

%% Parse inputs
if nargin<3,     help(mfilename); return, end
% if nargin<3,    n  = 2048; end
% if nargin<2,    bw = 2e3; end
% if nargin<1
%     p  = [  5,   1,   1,  10,    1,    1; ... 
%           921, 642, 453,   0, -528, -755; ... 
%           .1*ones(1,6)]; 
% end

f = linspace(-bw/2, bw/2, n);

%% Generate spectrum from free induction decay
N = size(p, 2);
A = p(1,:);  % amplitude
F = p(2,:);  % frequency
T = p(3,:);  % T2* (~= 1 / linewidth)

t   = (0:n-1)/bw;
fid = zeros(1, n);
for ii = 1:N
    if F(ii) > f(end)
        while F(ii) > f(end)
            F(ii) = (F(ii) - f(end)) + f(1); 
        end
        warning('genspec:aliasHi', ...
            'Peak at %g aliased to %g\n', p(2,ii), F(ii));
    elseif F(ii) < f(1)
        while F(ii) < f(1)
            F(ii) = (F(ii) - f(1)) + f(end);
        end
        warning('genspec:aliasLo', ...
            'Peak at %g aliased to %g\n', p(2,ii), F(ii));        
    end
    fid = fid + A(ii) * exp(-1i*2*pi * F(ii)*t) .* exp(-t/T(ii));
end
s = ifftdim(fid, 2);

% subplot 121;
% plot(t, real(fid), t, imag(fid));
% subplot 122;
% plot(f, real(s), f, imag(s), f, abs(s), 'k');
% set(gca, 'xdir', 'reverse')

