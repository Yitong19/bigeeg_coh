% Input directory
INDIR = 'C:\Users\yitong_z\Documents\BIGEEG\data\raw_data';
dir_real = 'C:\Users\yitong_z\Documents\BIGEEG\data\eeg_calcul\xcorr_permut';

% Parameters
origEventMarkers_resp = {'64532', '64533', '64534', '64535'};
EventMarkers_resp = {'Resp_RL_CCW', 'Resp_RL_CW', 'Resp_LR_CW', 'Resp_LR_CCW'};
Epoch_win = [-3, 7];                

% Reads all folders that are in INDIR
d = dir(dir_real);
isub = [d(:).isdir]; % returns logical vector if it's a folder
subjects = {d(isub).name}';
subjects(ismember(subjects, {'.', '..',})) = []; % Removes '.' and '..'

if isempty(subjects)    
    fprintf(sprintf('\n WARNING : NO DATA to process in %s\n',INDIR));
    return;
end

% Starts Brainsorm if it is not already running
if ~brainstorm('status')
    brainstorm;
end

% Set a Protocol name 
ProtocolName = 'BIGEEG_per_trial'; 

% Get the protocol index of an existing protocol 
iProtocol = bst_get('Protocol', ProtocolName);

% If protocol does not exist, create it
if isempty(iProtocol)
    gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);  % Create protocol without loading default anatomy
    iProtocol = bst_get('Protocol', ProtocolName);  % Re-fetch the protocol index after creation
end

% Loop through all subjects
for iSubj = 1:length(subjects)
    
    % Skip if subject already exists in the database
    if ~isempty(bst_get('Subject', subjects{iSubj}, 1))
        fprintf('Subject %s already exists in Brainstorm database - SKIP ----\n', subjects{iSubj});
        % continue;
    end

    % Prepare EEG file name
    fname = dir(fullfile(INDIR, subjects{iSubj}, 'eeg', strcat(subjects{iSubj}, '.bdf')));

    % Prepare BDF file + filepath
    RawFile = fullfile(INDIR, subjects{iSubj}, 'eeg', fname.name);

    % Get bad channels, ICA projectors, and bad ICA components files
    BadChannels = fullfile(INDIR, subjects{iSubj}, 'eeg', strcat(subjects{iSubj}, '_bad_channels.txt'));
    ICAprojectors = fullfile(INDIR, subjects{iSubj}, 'eeg', strcat(subjects{iSubj}, '_projectors.mat'));
    BadICAcomponents = fullfile(INDIR, subjects{iSubj}, 'eeg', strcat(subjects{iSubj}, '_bad_ica.txt'));

    % Start a new report
    sFiles = [];
    bst_report('Start', sFiles);

    % Add subject to Brainstorm database
    [~, ~] = db_add_subject(subjects{iSubj}, [], 1, 0);
    panel_protocols('UpdateTree'); % Update Protocol in GUI

    % Process: Create link to raw file
    sFiles = bst_process('CallProcess', 'process_import_data_raw', sFiles, [], ...
        'subjectname', subjects{iSubj}, ...
        'datafile', {RawFile, 'EEG-BDF'}, ...
        'channelreplace', 1, ...
        'channelalign', 1, ...
        'evtmode', 'value');

    % Rename events
    sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
        'src',   cell2mat(strcat(origEventMarkers_resp,',')), ...
        'dest',  cell2mat(strcat(EventMarkers_resp,',')));
   
    % Mark bad channels (if available)
    if exist(BadChannels, 'file')
        chan_sFiles = bst_process('CallProcess', 'process_channel_setbad', sFiles, [], ...
            'sensortypes', fileread(BadChannels));
    else
        fprintf('%s does not exist  ------- No global BAD channels\n', BadChannels);
    end

    % Add EEG positions
    sFiles_chan = bst_process('CallProcess', 'process_channel_addloc', sFiles, [], ...
        'channelfile', {'', ''}, ...
        'usedefault', 'ICBM152: BioSemi 64 10-10', ...
        'fixunits', 1, ...
        'vox2ras', 1);

    % Re-reference EEG
    sFiles = bst_process('CallProcess', 'process_eegref', sFiles, [], ...
        'eegref', 'AVERAGE', ...
        'sensortypes', 'EEG');

    % Process: High-pass:0.1Hz
    sFiles = bst_process('CallProcess', 'process_bandpass', sFiles, [], ...
        'sensortypes', 'EEG', ...
        'highpass',    0.1, ...
        'lowpass',     0, ...
        'tranband',    0, ...
        'attenuation', 'strict', ...  % 60dB
        'ver',         '2019', ...  % 2019
        'mirror',      0, ...
        'read_all',    1);

    % Handle ICA components
    if exist(BadICAcomponents,'file')
    
        % Load ICA projector file 
        projectors = load(ICAprojectors) ; 
    
        % Select ICA (picard) projector
        ica_proj = projectors.Projector(contains({projectors.Projector.Comment},'ICA_picard')); 
        
        % Select RE-Ref (AVERAGE) projector
        reref_proj = projectors.Projector(contains({projectors.Projector.Comment},'EEG reference: AVERAGE')); 
        
        % Reads the ICA components to reject 
        idx_to_reject = fileread(BadICAcomponents);
       
        % Get projector wich has 64 components (in case there are several tests with
        % various number of components we select the one that has 64
        for iIca=1:length(ica_proj)
            if size(ica_proj(iIca).Components,2)==64 ; ica_proj_64 = ica_proj(iIca); end
        end
    
        % Identify (flag) the components to remove
        ica_proj_64.CompMask = zeros(length(ica_proj_64.CompMask),1) ; 
        ica_proj_64.CompMask(str2num(idx_to_reject)) = 1 ;
    
        % Save projectors
        projectors.Projector(1) = reref_proj ; 
        projectors.Projector(2) = ica_proj_64 ; 
     
        % Get channel file 
        channelFile = bst_get('ChannelFileForStudy',  sFiles.FileName); 
        sChannel = in_bst_channel(channelFile);       
    
        % Append ICA projectors in the channel file 
        % sChannel.Projector = [sChannel.Projector, projectors.Projector] ; 
        sChannel.Projector = projectors.Projector;
        
        % Save channel file (will ICA projection apply automatically) 
        bst_save(file_fullpath(channelFile), sChannel, 'v7', 1);
      
    else 
        fprintf(strcat('YOU NEED TO IDENTIFY WICH COMPONENT TO REMOVE for subject : ',subjects{iSubj},'\n'));
        return ;       
    end

    % Interpolate bad electrodes
    sFiles = bst_process('CallProcess', 'process_eeg_interpbad', sFiles, [], ...
        'maxdist', 5, ...
        'sensortypes', 'EEG');

    % Process: Resample: 200Hz
    sFiles = bst_process('CallProcess', 'process_resample', sFiles, [], ...
        'freq',     200, ...
        'read_all', 1);

    % Loop through conditions response-locked
    for i_resp=1:length(EventMarkers_resp)
        
        % Import MEG/EEG: Events
        fprintf(EventMarkers_resp{i_resp})
    
        % Process: Import MEG/EEG: Events
        sFiles_epoched = bst_process('CallProcess', 'process_import_data_event', sFiles, [], ...
            'subjectname',   subjects{iSubj}, ...
            'condition',     '', ...
            'eventname',     EventMarkers_resp{i_resp}, ...
            'timewindow',    [], ...
            'epochtime',     Epoch_win, ...
            'split',         0, ...
            'createcond',    0, ...
            'ignoreshort',   1, ...
            'usectfcomp',    1, ...
            'usessp',        1, ...
            'freq',          [], ...
            'baseline',      [], ...
            'blsensortypes', 'MEG, EEG');

        % for iFile = 1:length(sFiles_epoched)
        %     % Get the current epoch name
        %     epoch_name = sFiles_epoched(iFile).Comment;
        % 
        %     % Modify the epoch name to create a new TF name
        %     new_name = strrep(epoch_name, '(#', '_');
        %     new_name = strrep(new_name, ')', '');
        % 
        %     % Rename the TF file using the new name
        %     sFiles_epoched(iFile) = bst_process('CallProcess', 'process_set_comment', sFiles_epoched(iFile), [], ...
        %         'tag', new_name, ...
        %         'isindex', 1);        
        % end
    end
end