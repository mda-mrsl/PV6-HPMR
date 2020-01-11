function p = cphantom(n)
%CPHANTOM makes complex Shepp-Logan phantom with non-linear phase
%   https://users.fmrib.ox.ac.uk/~mchiew/docs/PartialFourier.html
%
%   Usage: p = cphantom(n)
%
%       where n is the size of the phantom
%
%   See also PHANTOM
%
%   06/2019, Keith Michel

if ~nargin, help(mfilename), return; end

x  = (1:n); 
y  = (1:n)'; 
nn = round((n*[6, 6, 3, 9, 6])./16); 
im = 2*pi*(3*exp(-sqrt((x-nn(1)).^2+(y-nn(2)).^2)/nn(3))+...
            2*exp(-((x-nn(4)).^2)/nn(5)^2));
p  = phantom(n);
p  = p .* exp(1j*im);
