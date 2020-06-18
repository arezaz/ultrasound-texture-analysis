function [] = Func_StatisticalAnalysis(alpha)% Texture Analysis Code

% Statsitical Analysis Code
%% define read path to pooled folder
save_path = '.\Results';  % path for saving results and figures
pooled_path = [save_path, '\Pooled'];

%% loop over individual data and build pooled dataset
IndivFiles = dir([pooled_path, '/*_TextureAnalysisMetrics.csv']);  % read all csv files in pooled path
PooledResult = [];                          % pooled table init

disp(['-------- Stastical Analysis', ' --------']);

for k = 1:length(IndivFiles)
   
    read_from = [IndivFiles(k).folder, '/', IndivFiles(k).name];
    disp(['Reading: ', IndivFiles(k).name])
    IndivResult{k} = readtable(read_from);
    
    PooledResult = [PooledResult; IndivResult{k}];   
end

%% Save Pooled data
fn = 'TextureAnalysisMetrics_Pooled.csv';
writetable(PooledResult,[pooled_path, '\', fn]) %save pooled results to subject folder

%% Two-sample t-test btw

PooledisHealthyFalse = PooledResult(PooledResult.Subject_isHealthy == 0,:);
PooledisHealthyTrue = PooledResult(PooledResult.Subject_isHealthy == 1,:);
disp(['Pooling Reslts ...']);

%%
Features = PooledResult.Properties.VariableNames;
 
for fCount = 4:length(Features)
    i = fCount-3;
    Stats.Feature{i,:} = Features{fCount};
    [h,p,ci] = ttest2(PooledisHealthyTrue.(Features{fCount}),PooledisHealthyFalse.(Features{fCount}), 'Alpha',alpha);
    Stats.HypTestResult(i,:) = h;
    Stats.pval(i,:) = p;
end
disp(['Performing t-test ...']);

%% Save Stats Analysis
disp(['Saving Results in', save_path, ' ...']);

fn =  'STATS_TextureAnalysisMetrics_Pooled.csv';
writetable(struct2table(Stats),[pooled_path, '\', fn]) %save to subject folder
disp(['-------- Stastical Analysis Completed!', ' --------']);

