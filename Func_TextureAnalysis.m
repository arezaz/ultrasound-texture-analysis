function [] = Func_TextureAnalysis(filename,sub,ishealthy)% Texture Analysis Code
%% import dcm data
read_path = '.\DICOM_data\'; % data folder path
read_filename = filename; % file name DICOM format (.dcm)

save_path = '.\Results';  % path for saving results and figures

Subject = sub;  % Subject ID
isHealthy = ishealthy;  % healthy subject? true/false

read_dir = [read_path, read_filename];
save_dir = [save_path, '\Subject_', Subject];
mkdir(save_dir)

pooled_path = [save_path, '\Pooled'];
mkdir(pooled_path)
%% params for function input
Params.Subject = Subject;
Params.isHealthy = isHealthy;
Params.read_dir = read_dir;
Params.save_dir = save_dir;
Params.read_filename = read_filename;

%% Texture Analysis Funciton
Output = Texture(Params);

%% Save Metrics as Excel file
% Subject id and label
TextureMetrics.Subject_id = Output.Subject.ID;
TextureMetrics.Filename = Output.Subject.Filename;
TextureMetrics.Subject_isHealthy = Output.Subject.isHealthy;
% Convert to a single data structure
metricList = fieldnames(Output.Metrics);
for metric=1:numel(metricList)
    fn = fieldnames(Output.Metrics.(metricList{metric}));
    for field = 1:numel(fn)
        feature_name = [metricList{metric}, '_', fn{field}];
    TextureMetrics.(feature_name) = Output.Metrics.(metricList{metric}).(fn{field});
    end
end

% Save outputs as an Excel file
fn = [Output.Subject.ID, '_TextureAnalysisMetrics.csv'];
fn_structure = [save_dir,'/',Output.Subject.ID,'_TextureAnalysisMetrics.mat'];
save(fn_structure, 'TextureMetrics')
writetable(struct2table(TextureMetrics),[save_dir, '\', fn]) %save to subject folder
writetable(struct2table(TextureMetrics),[pooled_path, '\', fn]) %save to pooled folder

end