function n = nrmse(x, y)
%NRMSE computes normalized root mean square error
%
%   Usage: n = nscale(x, y)
%
%       where x is array with estimated or simulated data
%             y is array with ground truth data (only magnitude used)
%
%   See also NORM
%
%   06/2019, Keith Michel

%% Parse inputs
if nargin<2, help(mfilename); return; end

%% Compute NRMSE
y = abs(y);
n = norm(x(:) - y(:)) / norm(y(:));