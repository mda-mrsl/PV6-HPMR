function [x,n] = nscale(x)
%NSCALE scales data to range of [0, 1]
%
%   Usage: [y, n] = nscale(x)
%
%       where x is array with data to be rescaled
%             y is the normalized magnitude of x
%             n is a two-element vector containing the offset and scale
%               if x is real, y = (x - n(1)) / n(2)
%
%   See also SCAL, AUTOSC, UNSCALE
%
%   06/2019, Keith Michel

%% Parse inputs
if ~nargin,      help(mfilename); return; end
if isempty(x),   return; end
if ~isreal(x),   x = abs(x); end

%% Normalize
n(1) = min(x(:));
x    = x - min(x(:));
n(2) = max(x(:));
x    = x / max(x(:));
