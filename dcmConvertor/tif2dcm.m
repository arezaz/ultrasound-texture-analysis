function tif2dcm(ReadTifFilename,SaveDicomFilename)
dicomwrite(imread(ReadTifFilename), SaveDicomFilename)

end

