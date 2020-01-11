function x = unscale(x, n)
%UNSCALE scales data using specified offset and scale
%
%   Usage: x = unscale(y, n)
%
%       where y is array with data to be rescaled
%             x is the rescaled magnitude of y
%             n is a two-element vector containing the offset and scale
%               if y is real, x = n(2)*y + n(1)
%
%   See also SCAL, AUTOSC, NSCALE
%
%   08/2019, Keith Michel

%% Parse inputs
if ~nargin,      help(mfilename); return; end
if isempty(x),   return; end
if ~isreal(x),   x = abs(x); end

%% Rescale
x = x * n(2) + n(1);
