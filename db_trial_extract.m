% Input directory 
INDIR = 'C:\Users\yitong_z\Documents\BIGEEG\data\raw_data';

% Starts Brainsorm if it is not already running
if ~brainstorm('status')
    brainstorm;
end

% Read all folders that are in INDIR 
d = dir(INDIR); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Input files
sFiles = [];
SubjectNames = subjects;

% Define condition names
EventMarkers = {'RL_CCW','RL_CW','LR_CW','LR_CCW'};  % Add more if needed

% Define the base directory where you want to save the files
base_save_dir = 'C:\Users\yitong_z\Documents\BIGEEG\coherence_try';

% Initialize variables to store common data
Time = [];

% Loop through each subject
for si = 1:length(SubjectNames)
    subject_name = SubjectNames{si};
    fprintf('Processing subject: %s\n', subject_name);

    % Create the directory for this subject if it doesn't exist
    subject_dir = fullfile(base_save_dir, subject_name);
    if ~exist(subject_dir, 'dir')
        mkdir(subject_dir);
    end

    % Loop through each condition
    for ci = 1:length(EventMarkers)
        cond_name = EventMarkers{ci};
        var_name = sprintf('r_%s_%s', subject_name, cond_name);
        save_path = fullfile(subject_dir, [var_name '.mat']);
        
        % Check if the .mat file already exists to avoid reprocessing
        if exist(save_path, 'file')
            fprintf('File %s already exists - SKIP ----\n', save_path);
            continue;
        end

        fprintf('Processing condition: %s\n', cond_name);
        
        % Select time serie files for this condition and subject
        sFiles = bst_process('CallProcess', 'process_select_files_data', sFiles, [], ...
            'subjectname',   subject_name, ...
            'condition',     '', ...
            'tag',           cond_name, ...
            'includebad',    1, ...
            'includeintra',  0, ...
            'includecommon', 0, ...
            'outprocesstab', 'no'); 

        num_trials = length(sFiles);
        TS_cond = [];

        for trial = 1:num_trials
            TSraw = in_bst_data(sFiles(trial).FileName, 'F','Time');
            TS_cond(trial,:,:,:) = TSraw.F; 
        end
        
        data_struct.(var_name) = TS_cond;

        % Save the TF data for this condition directly to a .mat file
        save(save_path, '-struct', 'data_struct', '-v7.3');
        fprintf('Saved %s to %s\n', var_name, save_path);

        % Clear the variable to free up memory
        clear TS_cond data_struct;
    end
    Times = TSraw.Time; 
end

save(fullfile(base_save_dir, 'Times.mat'), 'Time');
