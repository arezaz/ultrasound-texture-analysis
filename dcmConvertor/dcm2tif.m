function dcm2tif(ReadDicomFilename,SaveTifFilename)
I = dicomread(ReadDicomFilename);
imwrite(I, SaveTifFilename);
end

