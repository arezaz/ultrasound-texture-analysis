%% dicom to tiff

ReadDicomFilename = 'HealthyFile-Original.dcm';
SaveTifFilename = 'HealthyFile-Original.tif';

dcm2tif(ReadDicomFilename,SaveTifFilename)
%% tiff to dicom

ReadTifFilename = 'HealthyFile-Original.tif';
SaveDicomFilename = 'HealthyFile-Reconstruct.dcm';

tif2dcm(ReadTifFilename,SaveDicomFilename)
%%