function rfOut = shiftSpsp(rfIn, grad, dz, t, gmr)

% Apply phase shift to spec-spat B1 pulse to achieve desired slice offset
%##########################################################################
%# rfOut = shiftSpsp(rfIn, grad, dz, t, gmr)
%#
%#  REQUIRED INPUTS:
%#      rfIn - Complex vector of spec-spat RF B1 waveform
%#        > Units are irrelevant, magnitude is passed unchanged
%#      grad - Vector of spec-spat gradient waveform        [mT/m]
%#      dz   - Slice position offset                        [mm]
%#
%#  OPTIONAL INPUT:
%#      t   - Timing resolution of RF and gradient waveforms   [us]
%#          Default = 8 us (default gradient resolution for Bruker)
%#      gmr - Gyromagnetic ratio                               [MHz/T]
%#          Default = 10.705 MHz/T (13C)
%#
%#  OUTPUTS:
%#      rfOut - Phase-modulated B1 waveform with specified slice offset
%#
%#  LITERATURE:
%#      Bernstein, King and Zhou. Handbook of MRI pulse sequences (2004). 
%#          p161
%#
%# Keith Michel
%# kamichel at mdanderson.org
%# 11/2018
%##########################################################################

if nargin<3, help(mfilename); return; end
funVersion = '20181112';

%% Parse inputs
if nargin<4, t = 8; end
if nargin<5, gmr = 10.705; end

%% Calculate excitation k-space and apply phase shift
% Note that by convention k integral is expressed backwards in time
k     = flip(gmr * t * cumtrapz(flip(grad(:))) * 1e-6);  % [cycles/mm]
rfOut = rfIn(:) .* exp(1i * 2*pi * k * dz);
