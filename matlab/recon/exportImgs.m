function exportImgs(imgs, outType, outPath, subject, method, acqp, reco)
% EXPORTIMGS writes image data to files
%
%   Usage: exportImgs([imgs, outType, outPath, subject, method, acqp, reco])
%
%       where imgs is the array of image data or string specifying ParaVision
%               scan directory to export
%               if complex, magnitude is used
%               if omitted, current directory is used
%             outType is one fo the following file types to write as output:
%               1, 'dicom' or 'dcm': DICOM (default)
%               2, 'tiff' or 'tif': TIFF
%               3, 'jpeg' or 'jpg': JPEG
%               4 or 'png': PNG
%               A 4-element logical vector can be used to specify multiple types
%             outPath is the directory into which output image files are written
%               if omitted, current directory is used
%             subject is the Paravision study subject header file structure
%             method is the Paravision scan method header file structure
%             acqp is the Paravision scan acqp header file structure
%             reco is the Paravision recon reco header file structure
%               readBrukerHeader() will be called for missing header structures
%
%   See also READBRUKERHEADER, VIEWBRUKERIMGS, RECOBRUKERKSPACE
%
%   10/2019, Keith Michel

funVersion = 'v20191119';

%% Parse inputs
if ~nargin
    if isdeployed, imgs = uigetdir('', 'Select Bruker Experiment Directory');
    else,          imgs = []; end
end
if isempty(imgs),    imgs = '.'; end
if isscalar(imgs),   imgs = num2str(imgs); end
if nargin<2,         outType = []; end
if isempty(outType), outType = 'dicom'; end
if nargin<3,         outPath = []; end
if isempty(outPath), outPath = '.'; end
if nargin<4,         subject = []; end
if nargin<5,         method = []; end
if nargin<6,         acqp = []; end
if nargin<7,         reco = []; end

if isdeployed
    fprintf('Entering %s (%s) at %s\n', mfilename, funVersion, datestr(now))
end

%% Fill out missing images and headers
if ischar(imgs)
    scanDir  = imgs;
    studyDir = fileparts(scanDir);
    [imgs, tmpReco, tmpAcqp] = viewBrukerImgs(scanDir, 1, false);
    if isempty(reco), reco = tmpReco; end
    if isempty(acqp), acqp = tmpAcqp; end
else
    scanDir  = '.';
    studyDir = '..';
end
if ~exist(fullfile(studyDir, 'subject'), 'file'), studyDir = '.'; end
if isempty(subject), subject = readBrukerHeader(fullfile(studyDir, 'subject')); end
if isempty(method),  method = readBrukerHeader(fullfile(scanDir, 'method')); end
if isempty(acqp),    acqp = readBrukerHeader(fullfile(scanDir, 'acqp')); end
if isempty(reco),    reco = readBrukerHeader(fullfile(scanDir, 'pdata', '1', 'reco')); end

%% Reshape images, set output type
assert(acqp.ACQ_dim > 1, 'exportImgs:1DScan', 'Cannot export 1D images');
[imgs,n] = nscale(imgs);
imgs     = reshape(imgs, size(imgs,1), size(imgs,2), []);
maxVal   = double(intmax('uint16'));
n(2)     = n(2) / maxVal;
imgs     = uint16(imgs * maxVal);
nImg     = size(imgs, 3);
baseName = sprintf('img%%0%dd', numel(num2str(nImg)));
bitDepth = 16;
if ischar(outType) || isscalar(outType)
    switch lower(outType)
        case {1, 'dicom', 'dcm'}
            out = [1 0 0 0];
        case {2, 'tiff', 'tif'}
            out = [0 1 0 0];
        case {3, 'jpeg', 'jpg'}
            out = [0 0 1 0];
        case {4, 'png'}
            out = [0 0 0 1];
    end
else
    out = outType;
end

%% Export images
outPath = fullfile(outPath, ...
    regexprep(subject.SUBJECT_name_string, '[\<\>\:\"\/\\|?*]', '_'), ...
    regexprep(subject.SUBJECT_study_name, '[\<\>\:\"\/\\|?*]', '_'), ...
    regexprep(acqp.ACQ_scan_name, '[\<\>\:\"\/\\|?*]', '_'));
if ~exist(outPath, 'dir'), mkdir(outPath); end

rot = reco.RECO_transposition;
rot = repmat(rot, 1, nImg/numel(rot));
xc  = method.PVM_SPackArrReadOffset;
xc  = repmat(xc, 1, nImg/numel(xc));
xw  = method.PVM_Fov(1);
yc  = method.PVM_SPackArrPhase1Offset;
yc  = repmat(yc, 1, nImg/numel(yc));
yw  = method.PVM_Fov(2);
zc  = method.PVM_SPackArrSliceOffset;
switch method.PVM_SpatDimEnum
    case '2D'
        sliceThick = method.PVM_SliceThick;
        sliceSepn  = method.PVM_SPackArrSliceDistance;
        zPosition  = zeros(1, sum(method.PVM_SPackArrNSlices));
        zw         = sliceSepn .* (method.PVM_SPackArrNSlices-1);
        idx        = 1;
        for ii = 1:numel(zw)
            nSlice = method.PVM_SPackArrNSlices(ii);
            zPosition(idx:idx+nSlice-1) = zc(ii) + ...
                linspace(-zw(ii)/2, zw(ii)/2, nSlice);
            idx = idx + nSlice;
        end
        sliceSepn  = repmat(sliceSepn, 1, nImg/numel(sliceSepn));
    case '3D'
        sliceThick = method.PVM_Fov(3) / method.PVM_Matrix(3);
        sliceSepn  = repmat(sliceThick, 1, nImg);
        zPosition  = zc + linspace(-method.PVM_Fov(3)/2, method.PVM_Fov(3)/2, ...
            method.PVM_Matrix(3));
end
zPosition = repmat(zPosition, 1, nImg/numel(zPosition));

% DICOM
if out(1)
    warning('off', 'images:dicom_add_attr:wrongAttribData');
    repTimes  = acqp.ACQ_repetition_time;
    repTimes  = repmat(repTimes, 1, nImg/numel(repTimes));
    echoTimes = acqp.ACQ_echo_time;
    echoTimes = repmat(echoTimes, 1, nImg/numel(echoTimes));
    excAngles = acqp.ACQ_flip_angle;
    excAngles = repmat(excAngles, 1, nImg/numel(excAngles));
    gradMat   = method.PVM_SPackArrGradOrient;
    gradMat   = repmat(gradMat, nImg/size(gradMat,1), 1, 1);
    switch lower(acqp.ACQ_patient_pos)
        case 'head_supine'
            ptPos = 'HFS';
        case 'head_prone'
            ptPos = 'HFP';
        case 'head_left'
            ptPos = 'HFL';
        case 'head_right'
            ptPos = 'HFR';
        case 'foot_supine'
            ptPos = 'FFS';
        case 'foot_prone'
            ptPos = 'FFP';
        case 'foot_left'
            ptPos = 'FFL';
        case 'foot_right'
            ptPos = 'FFR';
    end
    metaStruct = struct('Modality', 'MR', ...
        'Manufacturer', subject.ORIGIN, ...
        'InstitutionName', acqp.ACQ_institution, ...
        'ReferringPhysicianName', acqp.ACQ_operator, ...
        'SeriesDescription', acqp.ACQ_scan_name, ...
        'ManufacturerModelName', acqp.INSTRUM, ...
        'PatientName', subject.SUBJECT_name_string, ...
        'PatientID', subject.SUBJECT_id, ...
        'PatientSpeciesDescription', subject.SUBJECT_type, ...
        'PatientPosition', ptPos, ...
        'ScanningSequence', 'RM', ...
        'MRAcquisitionType', method.PVM_SpatDimEnum, ...
        'SequenceName', acqp.ACQ_method, ...
        'SliceThickness', sliceThick, ...
        'NumberOfAverages', acqp.NA * acqp.NAE, ...
        'ImagingFrequency', acqp.BF1, ...
        'ImagedNucleus', acqp.NUCLEUS, ...
        'NumberOfPhaseEncodingSteps', acqp.ACQ_size(2), ...
        'EchoTrainLength', acqp.ACQ_rare_factor, ...
        'PercentPhaseFieldOfView', 100 * acqp.ACQ_fov(2) / acqp.ACQ_fov(1), ...
        'PixelBandwidth', 2 * acqp.SW_h / acqp.ACQ_size(1), ...
        'SoftwareVersion', acqp.ACQ_sw_version, ...
        'ProtocolName', acqp.ACQ_scan_name, ...
        'StudyID', subject.SUBJECT_study_name, ...
        'ImagesInAcquisition', nImg, ...
        'RescaleType', 'US');

    for ii = 1:nImg
        curGradMat                         = squeeze(gradMat(ii,1:2,:)).';
        metaStruct.RepetitionTime          = repTimes(ii);
        metaStruct.EchoTime                = echoTimes(ii);
        metaStruct.SpacingBetweenSlices    = sliceSepn(ii);
        metaStruct.SliceLocation           = zPosition(ii);
        metaStruct.FlipAngle               = excAngles(ii);
        metaStruct.ImageOrientationPatient = curGradMat(:);
        metaStruct.RescaleIntercept        = n(1);
        metaStruct.RescaleSlope            = n(2);
        metaStruct.ImagePositionPatient    = ...
            [xc(ii)-xw/2, yc(ii)-yw/2, zPosition(ii)];
        switch rot(ii)
            case 0
                metaStruct.InPlanePhaseEncodingDirection = 'COL' ;
                metaStruct.PixelSpacing                  = method.PVM_SpatResol(1:2);
            case 1
                metaStruct.InPlanePhaseEncodingDirection = 'ROW' ;
                metaStruct.PixelSpacing                  = fliplr(method.PVM_SpatResol(1:2));
            case 2
                metaStruct.InPlanePhaseEncodingDirection = 'COL' ;
                metaStruct.PixelSpacing                  = method.PVM_SpatResol([1,3]);
            case 3
                metaStruct.InPlanePhaseEncodingDirection = 'BOTH' ;
                metaStruct.PixelSpacing                  = fliplr(method.PVM_SpatResol([2,3]));
        end
        fileName = sprintf([baseName, '.dcm'], ii);
        dicomwrite(imgs(:,:,ii), fullfile(outPath, fileName), ...
            metaStruct);
    end
end

% TIFF
if out(2)
    desc = sprintf('%s\n%s\n%s', subject.SUBJECT_name_string, ...
        subject.SUBJECT_study_name, acqp.ACQ_scan_name);
    for ii = 1:nImg
        fileName = sprintf([baseName, '.tif'], ii);
        imwrite(imgs(:,:,ii), fullfile(outPath, fileName), 'tif', ...
            'Description', desc);
    end
end

% JPEG
if out(3)
    desc = {subject.SUBJECT_name_string; subject.SUBJECT_study_name; ...
        acqp.ACQ_scan_name};
    tmp = nscale(double(imgs));
    for ii = 1:nImg
        fileName = sprintf([baseName, '.jpg'], ii);
        imwrite(tmp(:,:,ii), fullfile(outPath, fileName), 'jpg', ...
            'Comment', desc);
    end
end

% PNG
if out(4)
    desc = sprintf('%s\n%s\n%s', subject.SUBJECT_name_string, ...
        subject.SUBJECT_study_name, acqp.ACQ_scan_name);
    for ii = 1:nImg
        if rot(ii)
            xRes = method.PVM_Fov(1) / size(imgs,1) / 1e3;
            yRes = method.PVM_Fov(2) / size(imgs,2) / 1e3;
        else
            xRes = method.PVM_Fov(2) / size(imgs,1) / 1e3;
            yRes = method.PVM_Fov(1) / size(imgs,2) / 1e3;
        end
        fileName = sprintf([baseName, '.png'], ii);
        imwrite(imgs(:,:,ii), fullfile(outPath, fileName), 'png', ...
            'BitDepth', bitDepth, 'Description', desc, 'ResolutionUnit', 'meter', ...
            'Source', 'MRI', 'XResolution', xRes, 'YResolution', yRes);
    end
end

if isdeployed
    fprintf('Leaving %s at %s\n', mfilename, datestr(now))
end

end
