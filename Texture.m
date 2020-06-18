function [Out] = Texture(Params)

id = Params.Subject;
isHealthy = Params.isHealthy;
dir = Params.read_dir;
savedir = Params.save_dir;

% This funciton performs texture analysis for Ultrasound images
%   Input:
%       id (str): Subject ID 
%       dir (str): directory of the DICOM input image file
%       issave (bool): save the figures?
%       savedir (str): save images as .fig in this directory
%
%   Output:
%         - Image:
%             . RawData: Raw pixel values
%             . ROI: Region of interes pixel values
%             . ROI_detrend: ROI with background trend correction
%             . Mask: Mask created to trim down ROI
%             . ROI_Masked: ROI after subtracting Mask
%             . ROI_Masked_Pixelcount: Number of pixels in Masked ROI
%         - Metrics:
%             . MovingAvgFilter: Filtered with a moving average window
%                   .Stdev: Standard deviation filter
%                   .Range: Range filter 
%                   .Entropy: Entropy filter 
%             . FirstOrderStats: First-order surface analysis
%                   .Skew_Biased:
%                          .Values: For all unmasked pixels
%                          .Average: Averaged across unmasked pixels
%                   .Skew_unBiased: 
%                          .Values
%                          .Average
%                   .Kurtosis_Biased: 
%                          .Values
%                          .Average
%                   .Kurtosis_unBiased:
%                          .Values
%                          .Average
%                   .Entropy: 
%                          .Values
%                          .Average
%             . SecondOrderStats:
%                   .autoc: Autocorrelation
%                   .contr: Contrast
%                   .corrm: Correlation: matlab
%                   .corrp: Correlation
%                   .cprom: Cluster Prominence
%                   .cshad: Cluster Shade
%                   .dissi: Dissimilarity
%                   .energ: Energy: matlab
%                   .entro: Entropy
%                   .homom: Homogeneity: matlab
%                   .homop: Homogeneity
%                   .maxpr: Maximum probability
%                   .sosvh: Sum of sqaures: Variance
%                   .savgh: Sum average
%                   .svarh: Sum variance
%                   .senth: Sum entropy
%                   .dvarh: Difference variance
%                   .inf1h: Informaiton measure of correlation1 
%                   .inf2h: Informaiton measure of correlation2
%                   .homom: Inverse difference (INV) is homom
%                   .indnc: Inverse difference normalized (INN)
%                   .idmnc: Inverse difference moment normalized
%             . GaborFilter:
%                   .W(alpha)O(theta): For Wavelength alpha Direction theta
%                          .Values: For all unmasked pixels
%                          .Average: Averaged across unmasked pixels
%
% Alireza Rezazadeh
% rezaz003@umn.edu
% Spring 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['-------- Texture Analysis for Subject: ', id, ' --------']);
    disp(['Extracting DICOM data from:', dir]);
%%%%%%%%%%%% dicom extract info
info = dicominfo(dir);
    if info.ColorType == 'grayscale'
        X = im2double(dicomread(info)); % convert to double
    else
        X = rgb2gray(im2double(dicomread(info)));   % assert that image is in grayscale
    end 
%%%%%%%%%%%%

%%%%%%%%%%%% plot raw image
f_raw = figure('Name','Raw Image','NumberTitle','off');
figure(f_raw)
    imshow(X,[]);
    title('Raw Image')

%     if issave
%         disp(['Saving the figure to:', savedir])
%         saveas(figure(f_raw),[savedir,'/',id,'_Raw.fig']);
%     end

%%%%%%%%%%%% continue?
IsCnt;
%%%%%%%%%%%% save figure?
IsSaveFig(f_raw, id, '_Raw', savedir);
%%%%%%%%%%%%



%%%%%%%%%%%% crop ROI
f_crop = figure('Name','Crop ROI','NumberTitle','off');
figure(f_crop)
    hold on, title('Crop ROI')
    [X_crp, ~] = imcrop(X, []);
figure(f_crop)
    imshow(X_crp,[]);
    title('Cropped ROI')
    
%     if issave
%         disp(['Saving the figure to:', savedir])
%         saveas(figure(f_crop),[savedir,'/',id,'_Raw_ROI.fig']);
%     end
%%%%%%%%%%%%

%%%%%%%%%%%% continue?
IsCnt;
%%%%%%%%%%%% save figure?
IsSaveFig(f_raw, id, '_Raw_ROI', savedir);
%%%%%%%%%%%%

%%%%%%%%%%%% detrend ROI background
Tr =  detrend_2d(X_crp);
X_crp_d = X_crp - Tr;
%%%%%%%%%%%%

%%%%%%%%%%%% plot ROI surface
f_surf = figure('Name','ROI Surface Info','NumberTitle','off');
figure(f_surf)
%3d
    subplot(231)
        surf(X_crp)
        axis square, title('3D ROI')
    subplot(232)
        surf(Tr);
        axis square, title('3D Background Trend')
    subplot(233)
        surf(X_crp_d)
        axis square, title('3D ROI-Detrend')
%2d        
    subplot(234)
        imshow(X_crp, [])
        axis square, title('2D ROI')
    subplot(235)
        imshow(Tr, []);
        axis square, title('2D Background Trend')
    subplot(236)
        imshow(X_crp_d, [])
        axis square, title('2D ROI-Detrend')
        
%     if issave
%         disp(['Saving the figure to:', savedir])
%         saveas(figure(f_surf),[savedir,'/',id,'_ROI_SurfaceInfo.fig']);
%     end 
%%%%%%%%%%%%     

%%%%%%%%%%%% continue?
IsCnt;
%%%%%%%%%%%% save figure?
IsSaveFig(f_surf, id, '_ROI_SurfaceInfo', savedir);
%%%%%%%%%%%%

%%%%%%%%%%%% Mask for ROI
X_crp_d_mask = X_crp_d;

figure
iscnt = 1;
while iscnt == 1
    imshow(X_crp_d_mask, []);
        title('Crop Mask Area')
    h = drawfreehand; %draw something
    
    % crop out/in?
    cnt = questdlg('Marked Area:', ...
	'Which Area to Keep?', ...
	'Keep','Remove','Keep');
    % Handle response
    switch cnt
        case 'Keep'
            M = ~h.createMask();
            X_crp_d_mask(M) = NaN;
        case 'Remove'
            M = h.createMask();
            X_crp_d_mask(M) = NaN;
    end
    
    imshow(X_crp_d_mask, []);

    % continue crop?
    cnt = questdlg('Continue Crop?', ...
	'Crop', ...
	'Yes','Exit','Exit');
    % Handle response
    switch cnt
        case 'Yes'
            iscnt = 1;
        case 'Exit'
            iscnt = 0;
    end
end

Mask = ~isnan(X_crp_d_mask); % ultimate mask
%%%%%%%%%%%%     

%%%%%%%%%%%% plot masked ROI
f_mask = figure('Name','Mask ROI Area','NumberTitle','off');
figure(f_mask)
    subplot(311)
        imshow(X_crp_d, [])
        axis equal, title('ROI')
    subplot(312)
        imshow(Mask, [])
        axis equal, title('Mask')
    subplot(313)
        imshow(X_crp_d_mask, [])
        axis equal, title('Masked ROI')
        
%     if issave
%         disp(['Saving the figure to:', savedir])
%         saveas(figure(f_mask),[savedir,'/',id,'_ROI_Mask.fig']);
%     end
%%%%%%%%%%%%

%%%%%%%%%%%% continue?
IsCnt;
%%%%%%%%%%%% save figure?
IsSaveFig(f_mask, id, '_ROI_Mask', savedir);
%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Calculating Moving Average Window Filters ...');
%%%%%%%%%%%% Moving Average Window Filter-Based Metrics
kernel_size = 11;

X_crp_d_stdflt = stdfilt(X_crp_d, ones(kernel_size)); % Std Filter
X_crp_d_rngflt = rangefilt(X_crp_d, ones(kernel_size)); % Range Filter
X_crp_d_entflt = entropyfilt(X_crp_d, ones(kernel_size)); % Range Filter    
    
% plot filtered images
f_filt = figure('Name','Moving Window Filtered ROI','NumberTitle','off');
figure(f_filt)
    subplot(221)
        imshow(X_crp_d, [])
        axis equal, title('Original')
    subplot(222)
        imshow(X_crp_d_stdflt, [])
        axis equal, title('stdev filter')
    subplot(223)
        imshow(X_crp_d_rngflt, [])
        axis equal, title('range filter')
    subplot(224)
        imshow(X_crp_d_entflt, [])
        axis equal, title('entropy filter')    

        
%     if issave
%         disp(['Saving the figure to:', savedir])
%         saveas(figure(f_filt),[savedir,'/',id,'_ROI_FilteredImage.fig']);
%     end
%

%%%%%%%%%%%% continue?
IsCnt;
%%%%%%%%%%%% save figure?
IsSaveFig(f_filt, id, '_ROI_FilteredImage', savedir);
%%%%%%%%%%%%


% Metric Average for Masked ROI
metricList = {'Stdev', 'Range', 'Entropy'};
metricValueClass = {X_crp_d_stdflt(Mask), X_crp_d_rngflt(Mask), X_crp_d_entflt(Mask)}; %only keeping values for the unmasked area

numPix = nnz(~isnan(X_crp_d_mask)); % number of pixels in the masked ROI    

    for i = 1:length(metricList)
        MovingAvgFilter.(metricList{i}) = nansum(nansum(metricValueClass{i}))/numPix; %taking average by the coutns of pixels
    end

%%%%%%%%%%%%   

    disp('Calculating Gabor Filter ...');

%%%%%%%%%%%% Gabor Filter Texture Analysis
gaborArray = gabor([2 4 8 16],[0 45 90 135]); % wavelength and orientation

gaborMag = imgaborfilt(X_crp_d,gaborArray);

% plot filtered images
f_gabor = figure('Name','Gabor Filtered ROI','NumberTitle','off', 'Position', [55 55 1200 600]);
figure(f_gabor)
subplot(4,4,1);
for p = 1:16
    subplot(4,4,p)
    imshow(gaborMag(:,:,p),[]);
    theta = gaborArray(p).Orientation;
    lambda = gaborArray(p).Wavelength;
    title(sprintf('Orientation=%d, Wavelength=%d',theta,lambda));
end
        
%     if issave
%         disp(['Saving the figure to:', savedir])
%         saveas(figure(f_gabor),[savedir,'/',id,'_Gabor_FilteredImage.fig']);
%     end
%

%%%%%%%%%%%% continue?
IsCnt;
%%%%%%%%%%%% save figure?
IsSaveFig(f_gabor, id, '_Gabor_FilteredImage', savedir);
%%%%%%%%%%%%

% Metric Average for Masked ROI

numPix = nnz(~isnan(X_crp_d_mask)); % number of pixels in the masked ROI    

    for i = 1:16
        metricName = ['W' num2str(gaborArray(i).Wavelength) 'O' num2str(gaborArray(i).Orientation)];
        GaborFilter_temp = gaborMag(:,:,i);
        GaborFilter.(metricName) = nansum(nansum(GaborFilter_temp))/numPix; %only for unmasked area
    end
%%%%%%%%%%%%

    disp('Calculating First-Order Texture Metrics ...');

%%%%%%%%%%%% First-Order Statistical Analysis
    
X_crp_d_mask_unroll =  X_crp_d_mask(Mask); %only keepin the unmasked values

Skew_Biased= skewness(X_crp_d_mask_unroll(:));
Skew_unBiased = skewness(X_crp_d_mask_unroll(:),0);

Kurtosis_Biased =  kurtosis(X_crp_d_mask_unroll,1,'all');
Kurtosis_unBiased =  kurtosis(X_crp_d_mask_unroll,0,'all');

Entropy = entropy(X_crp_d_mask_unroll);    
    
% Metric Average for Masked ROI
metricList = {'Skew_Biased', 'Skew_unBiased', 'Kurtosis_Biased', 'Kurtosis_unBiased', 'Entropy'};
metricValueClass = {Skew_Biased, Skew_unBiased, Kurtosis_Biased, Kurtosis_unBiased, Entropy};

    for i = 1:length(metricList)
        FirstOrderStats.(metricList{i}) = metricValueClass{i};
    end

%%%%%%%%%%%%  

    disp('Calculating Second-Order Texture Metrics ...');

%%%%%%%%%%%% Second-Order Statistical Analysis
greyLevelNumber = 256;
glcms = graycomatrix(X_crp_d_mask,'NumLevels',greyLevelNumber); %calc glmatrics for masked roi

SecondOrderStats = GLCM_Features(glcms,0);
%%%%%%%%%%%%



%%%%%%%%%%%% Output Structured Data
disp('Generating Output Structure ...')
Out.Subject.ID = id;
Out.Subject.Filename = Params.read_filename;
Out.Subject.isHealthy = isHealthy;
% image data
Out.Image.RawData = X;
Out.Image.ROI = X_crp;
Out.Image.ROI_detrend = X_crp_d;
% filtered data
% Out.Image.Filtered.std = X_crp_d_stdflt;
% Out.Image.Filtered.range = X_crp_d_rngflt;
% Out.Image.Filtered.entropy = X_crp_d_entflt;
% masked image data
Out.Image.ROI_Masked = X_crp_d_mask;
Out.Image.Mask = Mask;
Out.Image.ROI_Masked_Pixelcount = numPix;
% metrics
Out.Metrics.MovingAvgFilter = MovingAvgFilter;      % moving average filter
Out.Metrics.FirstOrderStats = FirstOrderStats;      % first order analysis
Out.Metrics.SecondOrderStats = SecondOrderStats;    % second order analysis
Out.Metrics.GaborFilter = GaborFilter; % gabor filter

    %if issave % always save the output structure!
        disp(['Saving Output Data to:', savedir])
        filename = [savedir,'/',id,'_TextureAnalysis_Data.mat'];
        save(filename, 'Out')
    %end
    
    disp('-------- Texture Analysis Completed! --------');
end

