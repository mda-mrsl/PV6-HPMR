function studyString = readBrukerStudy( studyDir )
%READBRUKERSTUDY prints a study name and scan names for a ParaVision study folder
%
%   Usage: readBrukerStudy(studyDir)
%
%       where studyDir is a Bruker ParaVision study folder
%           if omitted, the working directory is used
%
%   See also READBRUKERHEADER, READBRUKERFID, RECOBRUKERKSPACE
%
%   09/2017, Chris Walker
%   06/2019, Keith Michel

if ~nargin, studyDir = '.'; end
subject = readBrukerHeader(fullfile(studyDir, 'subject'));
this    = sprintf('\nStudy %s\n\t', subject.SUBJECT_study_name);

subFolders = dir(studyDir);
scans      = zeros(numel(subFolders),1);
for ii =1:numel(subFolders)
    scans(ii) = str2double(subFolders(ii).('name'));
end
scans(isnan(scans)) = [];
scans = sort(scans);
for ii = 1:numel(scans)
    if exist(fullfile(studyDir,num2str(scans(ii)),'method'), 'file')
        m = readBrukerHeader(fullfile(studyDir,num2str(scans(ii)),'method'));
        a = readBrukerHeader(fullfile(studyDir,num2str(scans(ii)),'acqp'));
        s = a.ACQ_size;
        s(1) = s(1)/2;
        str = sprintf('%d x ', s);
        this = [this, ...
            sprintf('Scan # %d: %s, a (%s) x %d x %d %s Scan.\n\t', ...
            scans(ii), a.ACQ_scan_name, str(1:end-3), a.NSLICES, a.NR, m.Method)];
    else
        this = [this, ...
            sprintf('Scan No %d has no method header.\n\t',scans(ii))];
    end
end
this = sprintf('%s\n', this);

if ~nargout
    fprintf(this);
else
    studyString = this;
end

end

