%% This script is for preprocessing the behavioral data (removing pen lifts; resample)
% Input directory 
indir = 'C:\Users\yitong_z\Documents\BIGEEG\data\raw_data';
d = dir(indir);
isub = [d(:).isdir]; 
subject_names = {d(isub).name}';
subject_names(ismember(subject_names,{'.', '..'})) = []; 

flag = contains(subject_names,'bi');

freq_resample = 200;

out_dir = 'C:\Users\yitong_z\Documents\BIGEEG\data\beh_preprocess';

for isubj = 1:length(subject_names)

    out_path = fullfile(out_dir,[subject_names{isubj},'.mat']);

    if isfile(out_path)
        fprintf('Skipping %s: file already exists.\n', subject_names{isubj});
        continue
    end

    dir_beh = fullfile(indir,subject_names{isubj},'behavior');
    beh_list = dir(dir_beh);
    beh_blocks = {beh_list.name};
    beh_blocks(ismember(beh_blocks,{'.','..'})) = [];

    % Extract time stamps and put in the chronological order
    beh_order = cellfun(@(x)regexp(x,'\d{2}h\d{2}','match','once'),beh_blocks,'UniformOutput',false);

    times = datetime(beh_order,'InputFormat','HH''h''mm');
    [~, order] = sort(times);
    beh_blocks_sorted = beh_blocks(order);

    beh_clean = {};

    for iblock = 1:length(beh_blocks)

        block_path = fullfile(dir_beh,beh_blocks_sorted{iblock});

        txt = readtable(block_path);
        txt_clean = txt(:,{'participant','direction','rotation','trials_thisTrialN','stim_trial_started',...
                           'device_time','logged_time','time','x','y','z','pressure','orient_altitude','orient_azimuth'});
        txt_clean = renamevars(txt_clean,{'participant','trials_thisTrialN','stim_trial_started','device_time','pressure','orient_altitude','orient_azimuth'}, ...
                                         {'subj','tr_num','start','t','p','al','az'});

        trial_starts = [1;find(diff(txt_clean.tr_num)~=0)+1];
        trial_ends = [trial_starts(2:end)-1;height(txt_clean)];

        trial_tables = {};

        for ti = 1:numel(trial_starts)
            current_trial = txt_clean(trial_starts(ti):trial_ends(ti), :);

            if strcmp(current_trial.direction(1),'GD') && strcmp(current_trial.rotation(1),'antihoraire')
                current_trial.direction(:) = {'LR'};
                current_trial.rotation(:) = {'CCW'};
            elseif strcmp(current_trial.direction(1),'GD') && strcmp(current_trial.rotation(1),'horaire')
                current_trial.direction(:) = {'LR'};
                current_trial.rotation(:) = {'CW'};
            elseif strcmp(current_trial.direction(1),'DG') && strcmp(current_trial.rotation(1),'horaire')
                current_trial.direction(:) = {'RL'};
                current_trial.rotation(:) = {'CW'};
            elseif strcmp(current_trial.direction(1),'DG') && strcmp(current_trial.rotation(1),'antihoraire')
                current_trial.direction(:) = {'RL'};
                current_trial.rotation(:) = {'CCW'};
            else
                continue
            end

            % remove first pen lift
            if current_trial.p(1) == 0
                current_start = find(current_trial.p~=0,1,'first');
                current_trial = current_trial(current_start:end,:);
            end
            
            if size(current_trial,1) > 0
                % remove last if it exists
                if current_trial.p(end) == 0
                    current_trial_flip = flip(current_trial,1);
                    current_end = find(current_trial_flip.p~=0,1,'first');
                    current_trial_flip = current_trial_flip(current_end:end,:);
                    current_trial = flip(current_trial_flip,1);
                end
    
                % resample
                t_rel = double(current_trial.t);
                t_rel = t_rel - t_rel(1); 
                t_new_rel = (0:1/freq_resample:t_rel(end))';
    
                vars = {'x', 'y', 'z', 'p', 'al', 'az'};
    
                % Repeat the first row to match the new number of samples
                trial_resample = current_trial(ones(numel(t_new_rel), 1), :);
    
                % Replace time with the new relative timestamps
                trial_resample.t = t_new_rel;
    
                % Replace each numeric column with its interpolated values
                for vi = 1:numel(vars)
                    name = vars{vi};
                
                    trial_resample.(name) = interpolate_column( ...
                        t_rel, double(current_trial.(name)), t_new_rel);
                end
    
                trial_tables{end+1} = trial_resample; 
            else
                trial_tables{end+1} = current_trial;
            end
        end

        beh_clean{end+1} = trial_tables;
    end
    
    % save the data 
    save(out_path,'beh_clean');    
end


