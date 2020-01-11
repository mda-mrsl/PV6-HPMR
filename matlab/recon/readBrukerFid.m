function brukerFid = readBrukerFid(dirName)
%READBRUKERFID reads in raw MR data from Bruker scan directory
%
%   Usage: brukerFid = readBrukerFid(dirName)
%
%       where dirName is a Bruker scan directory
%               if omitted, the working directory is used
%             brukerFid is a vector containing the complex MR data
%               no reshaping or stripping of zeros is done
%
%   See also READBRUKERHEADER, READBRUKERSTUDY, RECOBRUKERKSPACE
%
%   06/2019, Keith Michel

%%
if ~nargin, dirName = '.'; end
dirName = num2str(dirName);

if exist(fullfile(dirName, 'fid'), 'file')
    f = fopen(fullfile(dirName, 'fid'), 'r');
elseif exist(fullfile(dirName, 'ser'), 'file')
    f = fopen(fullfile(dirName, 'ser'), 'r');
elseif exist(fullfile(dirName, 'rawdata.job0'), 'file')
    f = fopen(fullfile(dirName, 'rawdata.job0'), 'r');
    warning('readBrukerFid:jobAcq', ...
        'Specified directory uses acquisition jobs. Reading data for job 0.')
else
    error('readBrukerFid:fileNotFound', ...
        'Specified directory does not contain recognized raw data file')
end

%%
brukerFid = fread(f, inf, 'int32');
brukerFid = brukerFid(1:2:end) + 1i*brukerFid(2:2:end);
fclose(f);
