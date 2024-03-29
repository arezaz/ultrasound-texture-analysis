# Ultrasound Image Texture Analysis
Multiple texture analysis methods for ultrasound images in DICOM format. 
 
<img src="https://user-images.githubusercontent.com/45867828/150411986-91ada0f6-42ff-43f6-bbf2-67e6ccef3089.png" width="500">

### Citation
If you find this repository useful in your research, please consider citing the paper.

```
@misc{rezazadeh2022explainable,
      title={Explainable Ensemble Machine Learning for Breast Cancer Diagnosis based on Ultrasound Image Texture Features}, 
      author={Alireza Rezazadeh and Yasamin Jafarian and Ali Kord},
      year={2022},
      eprint={2201.07227},
      archivePrefix={arXiv},
      primaryClass={eess.IV}
}
```

## Instructions
### Texture Analysis for Ultrasound Images
  Input:  
  
       id (str): Subject ID  
       dir (str): directory of the DICOM input image file  
       issave (bool): save the figures and output data files?  
       savedir (str): save images as .fig in this directory  


   Output:  
   
         - Image:  
         
             . RawData: Raw pixel values  
             . ROI: Region of interes pixel values  
             . ROI_detrend: ROI with background trend correction  
             . Mask: Mask created to trim down ROI  
             . ROI_Masked: ROI after subtracting Mask  
             . ROI_Masked_Pixelcount: Number of pixels in Masked ROI 
             
         - Metrics:  
         
             . MovingAvgFilter: Filtered with a moving average window  
                   .Stdev: Standard deviation filter  
                   .Range: Range filter  
                   .Entropy: Entropy filter  
                   
             . FirstOrderStats: First-order surface analysis  
                   .Skew_Biased:  
                          .Values: For all unmasked pixels  
                          .Average: Averaged across unmasked pixels  
                   .Skew_unBiased:  
                          .Values  
                          .Average  
                   .Kurtosis_Biased:  
                          .Values  
                          .Average  
                   .Kurtosis_unBiased:  
                          .Values  
                          .Average  
                   .Entropy:  
                          .Values  
                          .Average  
                          
             . SecondOrderStats: Second-order surface analysis  
                   .autoc: Autocorrelation  
                   .contr: Contrast  
                   .corrm: Correlation: matlab  
                   .corrp: Correlation  
                   .cprom: Cluster Prominence  
                   .cshad: Cluster Shade  
                   .dissi: Dissimilarity  
                   .energ: Energy: matlab  
                   .entro: Entropy  
                   .homom: Homogeneity: matlab  
                   .homop: Homogeneity  
                   .maxpr: Maximum probability  
                   .sosvh: Sum of sqaures: Variance  
                   .savgh: Sum average  
                   .svarh: Sum variance  
                   .senth: Sum entropy  
                   .dvarh: Difference variance  
                   .inf1h: Informaiton measure of correlation1  
                   .inf2h: Informaiton measure of correlation2  
                   .homom: Inverse difference (INV) is homom  
                   .indnc: Inverse difference normalized (INN)  
                   .idmnc: Inverse difference moment normalized  
                   
             . GaborFilter: Gabor filter surface analysis  
                   .W(alpha)O(theta): For Wavelength alpha Direction theta  
                          .Values: For all unmasked pixels  
                          .Average: Averaged across unmasked pixels  

 Alireza Rezazadeh  
 rezaz003@umn.edu  
 Spring 2020  
