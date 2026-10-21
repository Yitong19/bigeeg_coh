% Add paths for toolboxes
tool_dir = 'C:\Users\yitong_z\Documents\MATLAB';
toolboxes = {'cartographie_motrice' 'cartographie_motrice/CM_extra' 'fieldtrip-20250107/external/fastica/' 'fieldtrip-20250107/external/mne/' 'c3d2OpenSim'};
for tool = toolboxes
    addpath(fullfile(tool_dir,tool{1}));
end
cd(fullfile(tool_dir,'fieldtrip-20250107'));
ft_defaults;

% Global file paths
dir_beh = 'C:\Users\yitong_z\Documents\BIGEEG\data\behavior_calcul\corrected_eeg';
dir_eeg = 'C:\Users\yitong_z\Documents\BIGEEG\data\eeg_calcul\reaction_times\corrected_times_no_outliers_beh';
dir_ts  = 'C:\Users\yitong_z\Documents\BIGEEG\coherence_try';
dir_out = 'C:\Users\yitong_z\Documents\BIGEEG\coherence_try\results';
if ~exist(dir_out, 'dir')
    mkdir(dir_out);
end

% Reads all folders that are in dir_beh 
d = dir(dir_beh);
isub = [d(:).isdir]; % returns logical vector if it's a folder
subject_names = {d(isub).name}';
subject_names(ismember(subject_names, {'.', '..'})) = []; % Removes '.' and '..'

flag = contains(subject_names,'bi');

% Define conditions 
Conds = {'RL_CCW', 'RL_CW', 'LR_CW', 'LR_CCW'};

% Load common structure data 
RowNames = load('C:\Users\yitong_z\Documents\BIGEEG\data\eeg_calcul\tf\RowNames.mat').RowNames;

% Common parameters
Fs = 200; % sampling frequency
Neeg = 64; % number of electrodes
window_seconds = 3; % window of analysis
L = round(Fs*window_seconds); % number of samples 

analysis_start_sample = 1;
analysis_idx = analysis_start_sample:(analysis_start_sample+L-1);

cmap = colormap(flipud(brewermap(1000,'RdBu')));

% Loop through all the subjects 
for i = 1:length(subject_names)

    % d = dir(dir_out); 
    % isub = [d(:).isdir]; % returns logical vector if is folder
    % subjects = {d(isub).name}';
    % subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..
    % 
    % if ismember(subject_names{i},subjects)
    %     fprintf('Subject %s already processed \n', subject_names{i});
    %     continue;
    % end 

    ts_path = fullfile(dir_ts,[subject_names{i}]);

    reac = fullfile(dir_eeg, [subject_names{i} '_reaction_time_beh_idx_cond.mat']);
    reaction_time = load(reac).beh_eeg_final;

    v_path = fullfile(dir_beh,subject_names{i},[subject_names{i},'_v.mat']);
    v = load(v_path).v_final;

    y_path = fullfile(dir_beh,subject_names{i},[subject_names{i},'_y.mat']);
    y = load(y_path).y_final;

    %% Regroup clean data by subject and by condition (cutting segments between stim and resp triggers)
    % Initialize frames for grouped clean data of tf and y pos (conditions ALL collapsed;
    % senses of rotation collapsed for LR and RL)
    ts_all = {};
    ts_RL = {};
    ts_LR = {};
    v_all = {};
    v_RL = {};
    v_LR = {};
    y_all = {};
    y_RL = {};
    y_LR = {};

    for c = 1:length(Conds)

        % Prepare the files for this condition
        reaction_time_cond = reaction_time.(Conds{c});

        v_cond = v.(Conds{c});
        v_cond_clean = cell2mat(v_cond(:,6:end));

        y_cond = y.(Conds{c});
        y_cond_clean = cell2mat(y_cond(:,6:end));

        ts = load(fullfile(ts_path,['r_' subject_names{i} '_' Conds{c}]));
        ts = ts.(['r_' subject_names{i} '_' Conds{c}]);

        % Reject the epochs
        eeg_cond = reaction_time_cond;
        rejectArray = zeros(size(eeg_cond,1),1);
        for t = 1:size(eeg_cond,1)
            if eeg_cond{t,4} == 0 
                rejectArray(t,1) = 2;
            elseif (isnan(eeg_cond{t,4}) ||  eeg_cond{t,5} == 0) && eeg_cond{t,4} ~= 0
                rejectArray(t,1) = 0;
            else
                rejectArray(t,1) = 1;
            end
        end

        if strcmp(subject_names{i}, 'FE37bi')
            rejectArray_stim = rejectArray;
            rejectArray_stim(rejectArray_stim==2) = 0;
            rejectArray(rejectArray==2) = [];           
        elseif strcmp(subject_names{i}, 'NO09')
            rejectArray_stim = rejectArray;
            rejectArray_stim(rejectArray_stim==2) = [];
            rejectArray(rejectArray==2) = 0;
        elseif strcmp(subject_names{i}, 'FE17') && strcmp(Conds{c}, 'RL_CW')
            rejectArray_stim = rejectArray;
            rejectArray_stim(end,1) = 0;
            rejectArray = rejectArray(1:end-1,:); 
            v_cond_clean = v_cond_clean(1:end-1,:);
            y_cond_clean = y_cond_clean(1:end-1,:);
        end

        reaction_time_valid = cell2mat(reaction_time_cond(logical(rejectArray),4));

        if sum(rejectArray)  == size(v_cond_clean,1) && ...
           size(v_cond_clean,1) == size(reaction_time_valid,1) && ...
           size(v_cond_clean,1) == size(y_cond_clean,1)

            ts_clean = ts(logical(rejectArray),:,:);
               
            [Ntr,Nch,Nt]=size(ts_clean);

            ts_whole = {};
            v_whole = {};
            y_whole = {};

            for tri = 1:Ntr
                time = floor((5-(reaction_time_valid(tri)))*200) - 2;
                ts_temp = {squeeze(ts_clean(tri,:,601:(601+time),:))};
                ts_whole = [ts_whole,ts_temp];
                v_temp = {squeeze(v_cond_clean(tri,201:(201+time)))};
                v_whole = [v_whole,v_temp];
                y_temp = {squeeze(y_cond_clean(tri,201:(201+time)))};
                y_whole = [y_whole,y_temp];
            end

            % Reorganize clean data into cells
            ts_all = [ts_all,ts_whole];
            v_all = [v_all,v_whole];
            y_all = [y_all,y_whole];

            if c < 3 % Conds (RL;RL;LR;LR)
                ts_RL = [ts_RL,ts_whole];
                v_RL = [v_RL,v_whole];
                y_RL = [y_RL,y_whole];
            else
                ts_LR = [ts_LR,ts_whole];
                v_LR = [v_LR,v_whole];
                y_LR = [y_LR,y_whole];
            end
        else
            fprintf('Need further verification of the matching for %s_%s \n', subject_names{i}, Conds{c});
        end
    end 

    Ntrials = numel(ts_all);
    
    ts_norm = cell(size(ts_all));
    v_norm = cell(size(v_all));
    y_norm = cell(size(y_all));
    
    for ti = 1:Ntrials
    
        ts = double(ts_all{ti});
        v  = double(reshape(v_all{ti}, 1, []));
        y  = double(reshape(y_all{ti}, 1, []));
    
        ts_norm{ti} = zscore(ts, 0, 2);
        v_norm{ti} = zscore(v);
        y_norm{ti} = zscore(y);
    end

    %% Cross correlation subject level (conditions are collapsed)
    total_len = cellfun(@(x) size(x,2),ts_norm);
    N_sep     = Ntrials - 1;
    N_total   = sum(total_len) + N_sep;

    ts_all_cat = zeros(Neeg,N_total);
    v_all_cat  = zeros(1,N_total);
    y_all_cat  = zeros(1,N_total);

    bad = false(1,N_total);
    trial_start = zeros(1,Ntrials);
    trial_end = zeros(1,Ntrials);
    separator_samples = [];

    current_sample = 1;

    for ti = 1:Ntrials
    
        eeg = double(ts_norm{ti});
        vel = double(reshape(v_norm{ti},1,[]));
        yel = double(reshape(y_norm{ti},1,[]));

        Nt = size(eeg,2);
    
        idx = current_sample:(current_sample+Nt-1);
    
        ts_all_cat(:,idx) = eeg;
        v_all_cat(idx) = vel;
        y_all_cat(idx) = yel;
    
        trial_start(ti) = idx(1);
        trial_end(ti) = idx(end);
    
        current_sample = idx(end)+1;
    
        if ti < Ntrials
    
            % Insert one artificial sample between trials
            ts_all_cat(:,current_sample) = 0;
            v_all_cat(current_sample) = 0;
            y_all_cat(current_sample) = 0;
            bad(current_sample) = true;
    
            separator_samples(end+1) = current_sample;
            current_sample = current_sample+1;
        end
    end

    data_cm = [];

    data_cm.label = RowNames(:);
    data_cm.fsample = Fs;
    
    data_cm.trial = {ts_all_cat};
    data_cm.time = {(0:N_total-1)/Fs};
    data_cm.sampleinfo = [1 N_total];
    data_cm.bad = bad;

    % Configurate coherence analysis
    fs = double(data_cm.fsample);

    CM = [];
    CM.quantum = round(2*fs);
    CM.tap = 3;
    CM.Find = 1:CM.quantum/Fs*20+1;
    CM.picksMEEG = 1:Neeg;
    CM.label = data_cm.label(CM.picksMEEG);

    % Use velocity as reference 
    CM.ref = struct();
    CM.ref(1).info   = 'velocity';
    CM.ref(1).chan   = double(v_all_cat(:)');
    CM.ref(1).filt   = [];
    CM.ref(1).norm   = 1;
    CM.ref(1).rectif = 0;

    CM.ref(2).info   = 'y_pos';
    CM.ref(2).chan   = double(y_all_cat(:)');
    CM.ref(2).filt   = [];
    CM.ref(2).norm   = 1;
    CM.ref(2).rectif = 0;

    CM = CM_coh_EEG_ref_one_pass(CM,data_cm,data_cm.bad);




   
    %% plot
    [max_v_ref,loc_v_ref] = findpeaks(squeeze(CM.Fyy(:,:,1)));  
    [~, order_v_ref] = sort(max_v_ref, 'descend');
    [max_y_ref,loc_y_ref] = findpeaks(squeeze(CM.Fyy(:,:,2))); 
    [~, order_y_ref] = sort(max_y_ref, 'descend');
    
    
    [max_v,loc_v] = findpeaks(squeeze(mean(CM.coh(1:64,:,1))));  
    [~, order_v] = sort(max_v, 'descend');
    
    [max_y,loc_y] = findpeaks(squeeze(mean(CM.coh(1:64,:,2)))); 
    [~, order_y] = sort(max_y, 'descend');
    
    fig = figure('Units','normalized', 'OuterPosition',[0 0 1 1], 'Color','w');

    t = tiledlayout(2,3, ...
        'TileSpacing','compact', ...
        'Padding','compact');
    
    nexttile(1);
    hData_v = plot(CM.f, squeeze(CM.Fyy(:,:,1))');
    set(hData_v, 'HandleVisibility', 'off');  
    hold on;
    
    if length(loc_v_ref) > 1
        f1_v = CM.f(loc_v_ref(order_v_ref(1)));
        f2_v = CM.f(loc_v_ref(order_v_ref(2)));
        
        h1_v = xline(f1_v, '--r', 'LineWidth', 1.5);
        h2_v = xline(f2_v, '--b', 'LineWidth', 1.5);
        
        legend([h1_v h2_v], ...
            {sprintf('%.2f Hz', f1_v), ...
             sprintf('%.2f Hz', f2_v)}, ...
            'Location', 'best');
    else
        f1_v = CM.f(loc_v_ref(order_v_ref(1)));
        
        h1_v = xline(f1_v, '--r', 'LineWidth', 1.5);
        
        legend(h1_v,sprintf('%.2f Hz', f1_v),'Location', 'best');
    end
    hold off;
    
    nexttile(4);
    hData_y = plot(CM.f, squeeze(CM.Fyy(:,:,2))');
    set(hData_y, 'Handlevisibility', 'off');  
    hold on;
    
    if length(loc_y_ref) > 1
        f1_y = CM.f(loc_y_ref(order_y_ref(1)));
        f2_y = CM.f(loc_y_ref(order_y_ref(2)));
        
        h1_y = xline(f1_y, '--r', 'LineWidth', 1.5);
        h2_y = xline(f2_y, '--b', 'LineWidth', 1.5);
        
        legend([h1_y h2_y], ...
            {sprintf('%.2f Hz', f1_y), ...
             sprintf('%.2f Hz', f2_y)}, ...
            'Location', 'best');
    else
        f1_y = CM.f(loc_y_ref(order_y_ref(1)));
        
        h1_y = xline(f1_y, '--r', 'LineWidth', 1.5);
        
        legend(h1_y,sprintf('%.2f Hz', f1_y),'Location', 'best');
    end
    hold off;
    
    
    
    nexttile(2);
    hData_vcoh = plot(CM.f, squeeze(CM.coh(1:64,:,1))');
    set(hData_vcoh, 'HandleVisibility', 'off');  
    hold on;
    
    f1_vcoh = CM.f(loc_v(order_v(1)));
    f2_vcoh = CM.f(loc_v(order_v(2)));
    
    h1_vcoh = xline(f1_vcoh, '--m', 'LineWidth', 1.5);
    h2_vcoh = xline(f2_vcoh, '--c', 'LineWidth', 1.5);
    
    legend([h1_vcoh h2_vcoh], ...
        {sprintf('%.2f Hz', f1_vcoh), ...
         sprintf('%.2f Hz', f2_vcoh)}, ...
        'Location', 'best');
    
    hold off;
    
    
    nexttile(5);
    hData_ycoh = plot(CM.f, squeeze(CM.coh(1:64,:,2))');
    set(hData_ycoh, 'Handlevisibility', 'off');  
    hold on;
    
    f1_ycoh = CM.f(loc_y(order_y(1)));
    f2_ycoh = CM.f(loc_y(order_y(2)));
    
    h1_ycoh = xline(f1_ycoh, '--m', 'LineWidth', 1.5);
    h2_ycoh = xline(f2_ycoh, '--c', 'LineWidth', 1.5);
    
    legend([h1_ycoh h2_ycoh], ...
        {sprintf('%.2f Hz', f1_ycoh), ...
         sprintf('%.2f Hz', f2_ycoh)}, ...
        'Location', 'best');
    
    hold off;
    
    cfg = [];
    cfg.xparam = 'time';
    cfg.layout = 'EEG1010.lay';
    cfg.figure = 'gca';
    cfg.comment = 'no';
    cfg.colormap = cmap;
    % cfg.interpolation = 'linear';
    % cfg.interplimits  = 'sensors';
    
    ave_v0 = [];
    ave_v0.avg     = squeeze(CM.coh(1:64,loc_v_ref(order_v_ref(1)),1));
    ave_v0.fsample = 1000;
    ave_v0.time    = 0;
    ave_v0.label   = data_cm.label(1:64);
    ave_v0.dimord  = 'chan_time';
    
    nexttile(3);
    ft_topoplotER(cfg, ave_v0);
    title(sprintf('%.2f Hz', CM.f(loc_v_ref(order_v_ref(1)))));
    
    % ave_v1 = [];
    % ave_v1.avg     = squeeze(CM.coh(1:64,loc_v(order_v_ref(2)),1));
    % ave_v1.fsample = 1000;
    % ave_v1.time    = 0;
    % ave_v1.label   = data_cm.label(1:64);
    % ave_v1.dimord  = 'chan_time';
    % 
    % nexttile(4);
    % ft_topoplotER(cfg, ave_v1);
    % title(sprintf('%.2f Hz', CM.f(loc_v(order_v_ref(2)))));
    
    ave_y0 = [];
    ave_y0.avg     = squeeze(CM.coh(1:64,loc_y_ref(order_y_ref(1)),2));
    ave_y0.fsample = 1000;
    ave_y0.time    = 0;
    ave_y0.label   = data_cm.label(1:64);
    ave_y0.dimord  = 'chan_time';
    
    nexttile(6);
    ft_topoplotER(cfg, ave_y0);
    title(sprintf('%.2f Hz', CM.f(loc_y_ref(order_y_ref(1)))));
    
    % ave_y1 = [];
    % ave_y1.avg     = squeeze(CM.coh(1:64,loc_y(order_y_ref(2)),2));
    % ave_y1.fsample = 1000;
    % ave_y1.time    = 0;
    % ave_y1.label   = data_cm.label(1:64);
    % ave_y1.dimord  = 'chan_time';
    % 
    % nexttile(8);
    % ft_topoplotER(cfg, ave_y1);
    % title(sprintf('%.2f Hz', CM.f(loc_y(order_y_ref(2)))));

    drawnow;

    output_file = fullfile(dir_out,sprintf('%s_coherence.png',subject_names{i}));

    saveas(fig, output_file);
    close(fig);

end
    % 

    % % Create the subject folder if it doesn't exist
    % subject_folder = fullfile(dir_out, subject_names{i});
    %     if ~exist(subject_folder, 'dir')
    %         mkdir(subject_folder);
    %     end
    % 
    % % Save the matrices
    % % mat_filename_all_unc = fullfile(subject_folder, sprintf('%s_all_unc.mat', subject_names{i}));
    % save(fullfile(subject_folder, sprintf('%s_all_unc.mat', subject)), ...
    %      '-fromstruct', 'struct', 'valcor_all', valcor_all);

    % %% Cross correlation subject + condition (LR)
    % tf_LR_cat = cat(2,tf_LR{:,:});
    % y_LR_cat = cat(2,y_LR{:,:});
    % y_LR_cat = y_LR_cat - mean(y_LR_cat);
    % valcor_LR = zeros(Nch,Nf,length(time_dur));
    % 
    % for chi=1:Nch
    %     for fi=1:Nf
    %         tf_signal = tf_LR_cat(chi,:,fi);
    %         tf_signal = tf_signal - mean(tf_signal);
    % 
    %         [valaux_LR,~] = xcorr(tf_signal, y_LR_cat, maxlag, 'coeff');
    %         valcor_LR(chi,fi,:) = valaux_LR;   
    %     end
    % end
    % 
    % % Save the matrices
    % % mat_filename_LR_unc = fullfile(subject_folder, sprintf('%s_LR_unc.mat', subject_names{i}));
    % % save(mat_filename_LR_unc, 'valcor_LR'); 
    % 
    % save(fullfile(subject_folder, sprintf('%s_LR_unc.mat', subject_names{i})), ...
    %      '-fromstruct', 'struct', 'valcor_LR', valcor_LR);
    % 
    % % Permutation analysis
    % % Dimension to permute: trials
    % valcor_LR_permuted = [];
    % 
    % for perm = 1:Nperm
    %     % Shuffle trial order
    %     permuted_indices = randperm(size(tf_LR,2));
    % 
    %     % Reinitialize permuted data storage
    %     y_whole_perm = [];
    % 
    %     % Concatenate trials with permuted order
    %     for tri = permuted_indices
    %         y_temp = y_LR{tri};
    %         y_whole_perm = cat(2, y_whole_perm, y_temp);
    %     end
    % 
    %     % Normalize permuted data
    %     y_norm_perm = y_whole_perm - mean(y_whole_perm);
    % 
    %     % Compute cross-correlation for permuted data
    %     for chi = 1:Nch
    %         for fi = 1:Nf
    %             tf_signal = tf_LR_cat(chi,:,fi);
    %             tf_signal = tf_signal - mean(tf_signal);
    % 
    %             [valaux_LR_perm,~] = xcorr(tf_signal, y_norm_perm, maxlag, 'coeff');
    %             valcor_LR_permuted(chi, fi, :, perm) = valaux_LR_perm;
    %         end
    %     end
    %     perm
    % end
    % 
    % % mat_filename_perm = fullfile(subject_folder, sprintf('%s_LR_perm.mat', subject_names{i}));
    % % save(mat_filename_perm, 'valcor_LR_permuted'); 
    % 
    % save(fullfile(subject_folder, sprintf('%s_LR_perm.mat', subject_names{i})), ...
    %      '-fromstruct', 'struct', 'valcor_LR_permuted', valcor_LR_permuted);
    % 
    % % Compute p-values
    % p_LR = zeros(Nch, Nf, length(time_dur)); % P-values
    % p_mask_LR = zeros(Nch, Nf, length(time_dur)); % Significance mask
    % valcor_LR_correct = zeros(Nch, Nf, length(time_dur)); % Corrected correlations
    % for chi = 1:Nch
    %     for fi = 1:Nf
    %         for t = 1:length(time_dur)
    %             % Compute mean and std of the permutation distribution
    %             mean_dist_LR = mean(valcor_LR_permuted(chi,fi,t,:),4);  
    %             std_dist_LR = std(valcor_LR_permuted(chi,fi,t,:), [], 4); 
    % 
    %             % Compute p-value using Gaussian approximation
    %             p_LR(chi, fi, t) = normcdf(valcor_LR(chi, fi, t), mean_dist_LR, std_dist_LR);
    % 
    %             % Apply significance mask (two-tailed)
    %             p_mask_LR(chi, fi, t) = (p_LR(chi, fi, t) < pval) + (p_LR(chi, fi, t) > (1 - pval));
    % 
    %             % Corrected correlation: keep only significant values
    %             valcor_LR_correct(chi, fi, t) = valcor_LR(chi, fi, t);
    %             if p_mask_LR(chi, fi, t) == 0
    %                 valcor_LR_correct(chi, fi, t) = NaN;
    %             end
    %         end
    %     end
    % end
    % 
    % % mat_filename_final = fullfile(subject_folder, sprintf('%s_LR_final.mat', subject_names{i}));
    % % save(mat_filename_final, 'valcor_LR_correct'); 
    % 
    % save(fullfile(subject_folder, sprintf('%s_LR_final.mat', subject_names{i})), ...
    %      '-fromstruct', 'struct', 'valcor_LR_correct', valcor_LR_correct);
    % 
    % %% Cross correlation subject + condition (RL)
    % tf_RL_cat = cat(2,tf_RL{:,:});
    % y_RL_cat = cat(2,y_RL{:,:});
    % y_RL_cat = y_RL_cat - mean(y_RL_cat);
    % valcor_RL = zeros(Nch,Nf,length(time_dur));
    % 
    % for chi=1:Nch
    %     for fi=1:Nf
    %         tf_signal = tf_RL_cat(chi,:,fi);
    %         tf_signal = tf_signal - mean(tf_signal);
    % 
    %         [valaux_RL,~] = xcorr(tf_signal, y_RL_cat, maxlag, 'coeff');
    %         valcor_RL(chi,fi,:) = valaux_RL;   
    %     end
    % end
    % 
    % % Save the matrices
    % % mat_filename_RL_unc = fullfile(subject_folder, sprintf('%s_RL_unc.mat', subject_names{i}));
    % % save(mat_filename_RL_unc, 'valcor_RL'); 
    % 
    % save(fullfile(subject_folder, sprintf('%s_RL_unc.mat', subject_names{i})), ...
    %      '-fromstruct', 'struct', 'valcor_RL_unc', valcor_RL_unc);
    % 
    % % Permutation analysis
    % % Dimension to permute: trials
    % valcor_RL_permuted = [];
    % 
    % for perm = 1:Nperm
    %     % Shuffle trial order
    %     permuted_indices = randperm(size(tf_RL,2));
    % 
    %     % Reinitialize permuted data storage
    %     y_whole_perm = [];
    % 
    %     % Concatenate trials with permuted order
    %     for tri = permuted_indices
    %         y_temp = y_RL{tri};
    %         y_whole_perm = cat(2, y_whole_perm, y_temp);
    %     end
    % 
    %     % Normalize permuted data
    %     y_norm_perm = y_whole_perm - mean(y_whole_perm);
    % 
    %     % Compute cross-correlation for permuted data
    %     for chi = 1:Nch
    %         for fi = 1:Nf
    %             tf_signal = tf_RL_cat(chi,:,fi);
    %             tf_signal = tf_signal - mean(tf_signal);
    % 
    %             [valaux_RL_perm,~] = xcorr(tf_signal, y_norm_perm, maxlag, 'coeff');
    %             valcor_RL_permuted(chi, fi, :, perm) = valaux_RL_perm;
    %         end
    %     end
    %     perm
    % end
    % 
    % % mat_filename_perm = fullfile(subject_folder, sprintf('%s_RL_perm.mat', subject_names{i}));
    % % save(mat_filename_perm, 'valcor_RL_permuted'); 
    % 
    % save(fullfile(subject_folder, sprintf('%s_RL_perm.mat', subject_names{i})), ...
    %      '-fromstruct', 'struct', 'valcor_RL_permuted', valcor_RL_permuted);
    % 
    % 
    % % Compute p-values
    % p_RL = zeros(Nch, Nf, length(time_dur)); % P-values
    % p_mask_RL = zeros(Nch, Nf, length(time_dur)); % Significance mask
    % valcor_RL_correct = zeros(Nch, Nf, length(time_dur)); % Corrected correlations
    % for chi = 1:Nch
    %     for fi = 1:Nf
    %         for t = 1:length(time_dur)
    %             % Compute mean and std of the permutation distribution
    %             mean_dist_RL = mean(valcor_RL_permuted(chi,fi,t,:),4);  
    %             std_dist_RL = std(valcor_RL_permuted(chi,fi,t,:), [], 4); 
    % 
    %             % Compute p-value using Gaussian approximation
    %             p_RL(chi, fi, t) = normcdf(valcor_RL(chi, fi, t), mean_dist_RL, std_dist_RL);
    % 
    %             % Apply significance mask (two-tailed)
    %             p_mask_RL(chi, fi, t) = (p_RL(chi, fi, t) < pval) + (p_RL(chi, fi, t) > (1 - pval));
    % 
    %             % Corrected correlation: keep only significant values
    %             valcor_RL_correct(chi, fi, t) = valcor_RL(chi, fi, t);
    %             if p_mask_RL(chi, fi, t) == 0
    %                 valcor_RL_correct(chi, fi, t) = NaN;
    %             end
    %         end
    %     end
    % end
    % 
    % % mat_filename_final = fullfile(subject_folder, sprintf('%s_RL_final.mat', subject_names{i}));
    % % save(mat_filename_final, 'valcor_RL_correct'); 
    % 
    % save(fullfile(subject_folder, sprintf('%s_RL_final.mat', subject_names{i})), ...
    %      '-fromstruct', 'struct', 'valcor_RL_correct', valcor_RL_correct);
% end
