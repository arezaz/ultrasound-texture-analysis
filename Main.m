%% Texture Analysis
filename = 'ser001img00512.dcm';
SubjectID = '015';
isHealthy = true;

Func_TextureAnalysis(filename,SubjectID,isHealthy)
%% Statistical Analysis
Alpha = 0.05; % t-test Significance Level
Func_StatisticalAnalysis(Alpha)

%%