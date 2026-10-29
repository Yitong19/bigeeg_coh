% ========================================================================
% This file is part of BIGEEG project
% 
% Free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This code is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%  
% Copyright (C) 2026 CNRS - Aix-Marseille Université
%
% ========================================================================
% This software was developed by
%       Yitong Zuo (Aix-Marseille Université)
% ------------------------------------------------------------------------
%
% DESCRIPTION : 
% Transform data so that they conform to the Brain Imaging Data Structure

% Input and output directories
indir    = 'C:\Users\yitong_z\Documents\BIGEEG\data\raw_data';
dir_bids = 'C:\Users\yitong_z\Documents\BIGEEG\data\BIDS';

% Originial and anonymized subject codes + meta data
sub_org = {'NO01bi','FE14','FE25','FE32bi','NO02bi','NO03','FE35bi','FE26bi','NO04',...
           'FE52bi','FE28bi','FE47bi','NO05bi','NO07','NO08','FE18','NO09','FE09','NO10',...
           'NO06','FE30bi','FE27','FE10bi','FE02','NO11','NO12bi','FE67bi','NO13bi','FE51bi',...
           'FE36bi','NO14','NO15bi','NO16bi','NO17bi','FE41bi','NO18','FE03','NO19','FE44bi',...
           'FE53bi','FE34bi','NO21','NO22bi','NO23bi','FE42bi','NO24bi','NO25bi','NO26bi',...
           'NO27bi','FE17','FE20','NO28','FE12','FE65','FE37bi','NO29','NO30','NO31','FE55',...
           'NO32','NO33','NO34'};
sub     = {'01','02','03','04','05','06','07','08','09','10','11','12','13',...
           '14','15','16','17','18','19','20','21','22','23','24','25','26','27','28',...
           '29','30','31','32','33','34','35','36','37','38','39','40','41','42','43',...
           '44','45','46','47','48','49','50','51','52','53','54','55','56','57','58',...
           '59','60','61','62'};
age     = [21 23 28 22 22 25 21 27 21 24 23 27 25 19 23 27 26 25 25 21 23 26 23 ...
           24 22 20 18 24 28 21 20 23 20 25 19 21 23 29 22 30 23 21 26 33 19 27 20 19 28 ...
           29 29 22 20 19 18 22 23 21 22 19 22 19];
sex     = {'f' 'f' 'f' 'm' 'm' 'f' 'm' 'f' 'f' 'f' 'f' 'f' 'm' 'm' 'm' 'm' 'm' 'f' ...
           'f' 'm' 'f' 'f' 'm' 'f' 'f' 'f' 'f' 'f' 'f' 'm' 'f' 'f' 'f' 'f' 'f' 'f' 'f' 'm' ...
           'f' 'm' 'm' 'f' 'f' 'f' 'm' 'm' 'f' 'f' 'f' 'm' 'm' 'f' 'm' 'f' 'f' 'f' 'm' 'f' ...
           'm' 'm' 'f' 'm'};
flag    = contains(sub_org,'bi');

% Behavioral json file
beh_var = {'subj' 'direction' 'rotation' 'block_num' 'tr_num' 'start' 't' 'logged_time' ...
           'time' 'x' 'y' 'z' 'p' 'al' 'az'};
descriptions = {
    'Subject identifier.'
    'Tracing direction.'
    'Tracing rotation.'
    'Block number.'
    'Trial number.'
    'Time stamp of trial commencement in seconds.'
    'Time stamp originating from the tablet (using its native timing reference), in milliseconds.'
    'Time stamp at which ioHub received the event, measured by PsychoPy clock, in seconds.'
    'Estimated time at which the event occurred in seconds.'
    'Pen tip position in the horizontal dimension.'
    'Pen tip position in the vertical dimension.'
    'Distance between the pen tip and the screen.'
    'Pressure exerted by the pen tip.'
    'Angle of tilt-altitute of the pen tip.'
    'Angle of tilt-alzimuth of the pen tip.'
};

jsonData = struct();

for i = 1:numel(beh_var)
    jsonData.(beh_var{i}).Description = descriptions{i};
end

jsonText = jsonencode(jsonData,'PrettyPrint',true);

% General information that applies to all
general= [];
general.InstitutionName                 = 'Aix-Marseille Université';
general.InstitutionalDepartmentName     = 'Centre de Recherche en Psychologie et Neurosciences, UMR 7077';
general.InstitutionAddress              = '3 place Victor Hugo, 13003 Marseille';
general.bidsroot = dir_bids;

% Loop through all the subjects
for subindx=1:numel(sub)

    %% EEG
    % Original file
    bdf_file = fullfile(indir,sub_org{subindx},'eeg',[sub_org{subindx},'.bdf']);
    
    % Save in BIDS form
    cfg_eeg = general;
    cfg_eeg.method    = 'convert';
    cfg_eeg.suffix    = 'eeg';
    cfg_eeg.dataset   = bdf_file;
    cfg_eeg.sub       = sub{subindx};
    cfg_eeg.participants.age = age(subindx);
    cfg_eeg.participants.sex = sex{subindx};
    cfg_eeg.participants.grp = flag(subindx);

    cfg_eeg.TaskName        = 'looptracing';
    cfg_eeg.TaskDescription = 'Subjects were tracing continuous loops (left-to-right/right-to-left; clockwise/counterclockwise)';
    
    cfg_eeg.eeg.PowerLineFrequency = 50; 
    cfg_eeg.eeg.SamplingFrequency  = 512; 
    cfg_eeg.eeg.EEGReference       = 'Average'; 
    
    data2bids(cfg_eeg);

    %% Behavior
    % Hard to use the FieldTrip function, save data as tsv files and the
    % json files created manually
    % Original folder
    dir_beh_out = fullfile(dir_bids,['sub-',sub{subindx}],'beh');
    if ~exist(dir_beh_out,'dir')
        mkdir(dir_beh_out);
    end
    dir_beh = fullfile(indir,sub_org{subindx},'behavior');
    beh_list = dir(dir_beh);
    beh_blocks = {beh_list.name};
    beh_blocks(ismember(beh_blocks,{'.','..'})) = [];
    
    % Extract time stamps and put in the chronological order
    beh_order = cellfun(@(x)regexp(x,'\d{2}h\d{2}','match','once'),beh_blocks,'UniformOutput',false);
    
    times = datetime(beh_order,'InputFormat','HH''h''mm');
    [~, order] = sort(times);
    beh_blocks_sorted = beh_blocks(order);
    
    % Initialize cell and store cleaned behavioral data 
    beh_clean = table();

    for iblock = 1:length(beh_blocks)

        block_num = iblock;

        block_path = fullfile(dir_beh,beh_blocks_sorted{iblock});

        txt = readtable(block_path);
        txt_clean = txt(:,{'participant','direction','rotation','trials_thisTrialN','stim_trial_started',...
                           'device_time','logged_time','time','x','y','z','pressure','orient_altitude','orient_azimuth'});
        txt_clean = renamevars(txt_clean,{'participant','trials_thisTrialN','stim_trial_started','device_time','pressure','orient_altitude','orient_azimuth'}, ...
                                         {'subj','tr_num','start','t','p','al','az'});

        trial_starts = [1;find(diff(txt_clean.tr_num)~=0)+1];
        trial_ends = [trial_starts(2:end)-1;height(txt_clean)];

        trial_tables = table();

        tri_num = 0;
        for ti = 1:numel(trial_starts)
            current_trial = txt_clean(trial_starts(ti):trial_ends(ti), :);

            % Rename condition in English
            if strcmp(current_trial.direction(1),'GD') && strcmp(current_trial.rotation(1),'antihoraire')
                current_trial.direction(:) = {'LR'};
                current_trial.rotation(:) = {'CCW'};
                tri_num = tri_num+1;
            elseif strcmp(current_trial.direction(1),'GD') && strcmp(current_trial.rotation(1),'horaire')
                current_trial.direction(:) = {'LR'};
                current_trial.rotation(:) = {'CW'};
                tri_num = tri_num+1;
            elseif strcmp(current_trial.direction(1),'DG') && strcmp(current_trial.rotation(1),'horaire')
                current_trial.direction(:) = {'RL'};
                current_trial.rotation(:) = {'CW'};
                tri_num = tri_num+1;
            elseif strcmp(current_trial.direction(1),'DG') && strcmp(current_trial.rotation(1),'antihoraire')
                current_trial.direction(:) = {'RL'};
                current_trial.rotation(:) = {'CCW'};
                tri_num = tri_num+1;
            else
                continue
            end

            current_trial(:,{'subj'}) = repmat(cellstr(sub{subindx}),size(current_trial,1),1);
            current_trial(:,{'tr_num'}) = array2table(repmat(tri_num,size(current_trial,1),1));
            current_trial.block_num = repmat(block_num,size(current_trial,1),1);
            current_trial = movevars(current_trial,"block_num",After="rotation");
            trial_tables = vertcat(trial_tables,current_trial);
        end
        beh_clean = vertcat(beh_clean,trial_tables);
    end

    beh_clean_path = fullfile(dir_beh_out,['sub-',sub{subindx},'_task-looptracing_beh.tsv']);
    writetable(beh_clean,beh_clean_path,'filetype','text','delimiter',',');

    beh_json_path = fullfile(dir_beh_out,['sub-',sub{subindx},'_task-looptracing_beh.json']);
    fid = fopen(beh_json_path,'w');
    fprintf(fid, '%s\n', jsonText);
    fclose(fid);
end

