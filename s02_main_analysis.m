
exp_dir = 'C:\Users\Antonella Iannotta\OneDrive - Université Libre de Bruxelles\Documents\data_DCD';
toolboxes = {'cartographie_motrice' 'cartographie_motrice/CM_extra' 'fieldtrip-20190205/external/fastica/' 'fieldtrip-20190205/external/mne/' 'c3d2OpenSim'};
for tool = toolboxes
    addpath(fullfile(exp_dir,'toolboxes',tool{1}));
end
cd(fullfile(exp_dir,'toolboxes','fieldtrip-20190205'));
ft_defaults;

exp_dir = 'C:\Users\Antonella Iannotta\OneDrive - Université Libre de Bruxelles\Documents\data_DCD';
cd(fullfile(exp_dir,'analysis'));
data_dir = fullfile(exp_dir,'eeg');
force_dir = fullfile(exp_dir,'force');
matfold = fullfile(exp_dir,'matfiles');
picsfold = fullfile(exp_dir,'pics');

subjects = [];
conds = {'EOhard' 'EOfoam' 'EChard' 'ECfoam'};
for n_sub = 1:52
    subjects{n_sub} = 'DCD000';
    this_num = num2str(n_sub);
    subjects{n_sub}(end-length(this_num)+1:end) = this_num;
end
% subjects([13 16 17 20:21 23 30 33 35 43]) = [];
subjects([13 16 17 20:21 23 30 33:52]) = []; % partcipants for CKC article adults vs kids (25 kids) - load this when needed


subjects_kids = subjects;

for n_sub = 1:62
    subjects{n_sub} = 'bal000';
    this_num = num2str(n_sub);
    subjects{n_sub}(end-length(this_num)+1:end) = this_num;
end
% subjects = subjects([1:5 7:13 15:20 25 27 28 31 36 37 39 42 45 46 60 62]);
subjects = subjects([1:5 7:13 16:20 25 27 28 31 36 37 39 42 45 46 60 62]); % removed 15
% subjects = subjects([26:35 37:39]);
subjects_adults = subjects;


subjects = cat(2,subjects_kids,subjects_adults);




SDcop = [];
Vcop = [];
 % for n_sub =  1
for n_sub =  26:length(subjects)
    
    for n_cond = 1
        mat_outfile = fullfile(matfold,[subjects{n_sub} '_' conds{n_cond} '_coherence.mat']);
        load(mat_outfile,'sway')
        % sway = sway(1);
        SDcop(n_sub,n_cond,:) = mean(cat(1,sway.SDcop_filt),1);
        Vcop(n_sub,n_cond,:) = mean(cat(1,sway.Vcop_filt),1);
    end
end  
   
 Vis_acc = [];
for n_sub = 1:length(subjects)
    vis_fold = fullfile(exp_dir,'vision');
    vis_file = fullfile(vis_fold,[subjects{n_sub} '_1.mat']);
    if exist(vis_file)
        load(vis_file)
        Vis_acc(n_sub) = result.Tfinal;
    else
        Vis_acc(n_sub) = NaN;
    end
end







%% Step 2 Coherence analysis

perm_stat = 0;
% load nave_corrected nave
% Max_ave = nave;
for n_sub = 1:length(subjects)

    subfold = fullfile(data_dir,subjects{n_sub});
    for n_cond = 1:length(conds)
      

        if strcmp(subjects{n_sub}(1:3),'bal')
            filename = [subjects{n_sub} '_part*_' conds{n_cond} '_ICA_new.mat'];
        else
            filename = [subjects{n_sub} '_part*_' conds{n_cond} '_ICA.mat'];
        end
        files = dir(fullfile(subfold,filename));


        
        if length(files) < 2
            disp(['Missing file for ' subjects{n_sub} '   ' conds{n_cond}])
        end
        
        Ndf = [];
        NdfSL = [];
        for n_file = 1:length(files)
            infile = fullfile(subfold,files(1).name);
            matfile = fullfile(subfold,files(n_file).name);
            this_data = load(matfile,'data','force');
            data = this_data.data;
            if strcmp(infile(end-27:end),'bal025_part11_EChard_ICA.mat')
                data.bad(:) = 0;
            end
            % if isfield(this_data,'force')
            force = this_data.force;
            Fs = 1000;
            Neeg = 61;

            deb = 10*Fs; fin = length(force)-10*Fs;
            % deb = 5*Fs; fin = 55*Fs;
            [tmp,filename] = fileparts(infile);
            if strcmp(filename,'DCD026_part1_EOhard_ICA')
                deb = 15*Fs;
            end
            cop = force([5 4],:).*[-1 1]'./force(3,:);
            cop = cop - nanmean(cop(:,deb:fin),2);
            cop(isnan(cop)) = 0;
            cop(:,1:deb) = cop(:,deb*ones(1,deb));
            cop(:,fin:end) = cop(:,fin*ones(1,length(cop)-fin+1));
            data.bad(1:deb) = 1;
            data.bad(fin:end) = 1;
            fig = figure('Menubar','none','Visible','Off');
            sway = CM_COP_freq_analysis(cop(:,deb:fin),Fs);
            print('-dpng',fullfile(picsfold,[subjects{n_sub} '_' conds{n_cond} '_file' num2str(n_file)]));
            close(fig)


            if n_file == 1
                all_sway = sway;
            else
                all_sway(n_file) = sway;
            end
                 
            if  0
                data.time{1}(1:200) = [];
                data.trial{1}(:,1:200) = [];
                data.sampleinfo(end) = data.sampleinfo(end)-200;
                data.bad(1:200) = [];
            end
            
            % degrees of freedom
            try
                [U,S,V] = svd(data.trial{1}(1:Neeg,find(~data.bad)),'econ');
                Ndf(n_file) = sum(diag(S)>S(1)/100);
            end

            surflap = 0;
            if surflap
                 cfg                 = [];
                 cfg.method          = 'finite';
                 cfg.elecfile        = 'standard_1005.elc';
                 cfg.elec            = ft_read_sens(cfg.elecfile);
                 cfg.trials          = 'all';
                 % cfg.neighbours      = neighbours;
                 cfg.conductivity    = 0.33;
                 dataSL               = ft_scalpcurrentdensity(cfg, data);
                 data.trial{1}(1:Neeg,:) = dataSL.trial{1}(1:Neeg,:);
            end
            
            % initialise CM
            surrog = [];
            CM = [];
            CM.win = 'boxcar';
            CM.quantum = 2*Fs;
            CM.tap = 3;
            CM.Find = 1:CM.quantum/Fs*10+1;
            
            
            F_h_x = cosine_filter(size(force,2),{'high' 'low'},[0.2 10]/Fs*2,[0.1 2]/Fs*2);
            cop_filt = real(ifft(fft(cop,[],2).*F_h_x,[],2));
            vcop = [0 cop_filt(2,3:end)-cop_filt(2,1:end-2) 0]*Fs/2;
            

            ref = [];
            ref.info = 'vcop'; % more afferent CKC
            ref.chan = vcop;
            ref.filt = [];
            ref.norm = 0;
            ref(2) = ref(1);
            ref(2).info = 'abs(vcop)';
            ref(2).chan = abs(vcop);
            ref(3) = ref(1); % more efferent CKC
            ref(3).info = 'abs(rcop)';
            ref(3).chan = abs(cop_filt(2,:));
            
            pct_vcop = prctile(vcop,[25 50 75]);
            d = diff(pct_vcop([1 3]));
            bad_cop = abs(vcop-pct_vcop(2)) > 5*d;
            bad_cop = conv(bad_cop,ones(1,Fs),'same') ~= 0;
            data.bad(find(bad_cop)) = 1;

            
            CM.infile = infile;
            CM.ref = ref;


            % surrogate-data-based stat
            if 1 % ~perm_stat & isfield(CM,'maxave')
                CM.surrog.Nsim = 1000;
                CM.surrog.Find = 1:21;            % 9 à 17 si analyse bande 4-8 (laisser 2 pour bande beta)
                % CM.surrog.sens = sensors.picksMEG;
            end
            
            CM.picksMEEG = 1:Neeg;
            CM.label = data.label(CM.picksMEEG) ;

            if exist('Max_ave')
                if n_file == 1
                    CM.maxave = Max_ave(n_sub,n_cond);
                else
                    CM.maxave = Max_ave(n_sub,n_cond)-CMall.nave;
                end
            end

            % if CM.maxave > 0
            % Truncate EEG signal to match force vector length

                CM = CM_coh_EEG_ref_one_pass(CM,data,data.bad);
    
                if n_file == 1 
                    CMall = CM;
                else
                    if CMall.nave == 0
                        CMall = CM;
                    else
                        if CM.nave > 0
                            CMall = CM_combine_results(CMall,CM);
                        end
                    end
                end
            % end
        end

        sway = all_sway;
        surrog = [];
        if isfield(CMall,'surrog')
            surrog = CMall.surrog;
            CMall.surrog = [];
        end
        
        mat_outfile = fullfile(matfold,[subjects{n_sub} '_' conds{n_cond} '_coherence.mat']);
        save(mat_outfile,'CMall','surrog','sway')
    end
end


% 
% %%%
% mat_outfile = fullfile(matfold,[subjects{n_sub} '_' conds{n_cond} '_coherence.mat']);
% r0 = load(mat_outfile,'CMall','surrog','sway')
% mat_outfile = fullfile(matfold,[subjects{n_sub} '_' conds{n_cond} '_coherence_new.mat']);
% r1 = load(mat_outfile,'CMall','surrog','sway')
% 
% n_ref = 1;
% figure
% plot(r0.CMall.f, r0.CMall.coh(1:61, :, n_ref)', 'k')
% figure
% n_ref = 1;
% figure
% plot(r0.CMall.f, r1.CMall.coh(1:61, :, n_ref)', 'k')
% %%%



Varcop = [];
for n_sub = 1
    for n_cond = 1:3
        mat_outfile = fullfile(matfold,subjects{n_sub},[subjects{n_sub} '_' conds{n_cond} '_coherence.mat']);
        load(mat_outfile)

        Varcop(:,n_cond,n_sub) = mean(CMall.Varcop,2);
    end
end

[h,p,tmp,stat] = ttest(permute(Varcop(:,[1 1 2],:)-Varcop(:,[2 3 3],:),[3 1 2]))






F = [3 5];
Find = round(F(1)/CMall.Fs*CMall.quantum+1):round(F(end)/CMall.Fs*CMall.quantum+1);
[M,n_elec] = max(mean(CMall.coh(1:Neeg,Find,n_ref),2));

subplot(2,2,1)
plot(CMall.f,(CMall.coh(n_elec,:,n_ref)))
xlabel('Frequency [Hz]'); ylabel('Coherence')

subplot(2,2,3)
plot(CMall.f,phase(CMall.Fxy(n_elec,:,n_ref)),'o')
xlabel('Frequency [Hz]'); ylabel('Phase')

subplot(2,2,2)
plot((-CMall.quantum/2+1:CMall.quantum/2)/Fs*1000,CMall.fxy(n_elec,:,n_ref))
xlabel('Time [ms]'); title('cross-correlogram')




plot(CMall.fxy(n_elec,:,n_ref))


plot(CMderivall.f,CMderivall.coh(1:Nderiv,:,n_ref)')


n_ref = 1;
n_f = find(CMderivall.f <= 5 & CMderivall.f >= 3);
to_plot_deriv = mean(CMderivall.coh(1:Nderiv,n_f,n_ref),2);
to_plot = zeros(Neeg,1);
for n_elec = 1:Neeg
    in_deriv = find(any(all_derivs' == n_elec));
    to_plot(n_elec) = mean(to_plot_deriv(in_deriv));
end


if 0
    cfg = [];
    cfg.xparam = 'time';
    % cfg.zlim = [0 0.16]
    cfg.layout = 'EEG1010.lay';
    ave = [];
    ave.avg = to_plot;
    cfg.zlim = [0 max(ave.avg)];
    cfg.comment = 'no';
    ave.fsample = 1000;
    ave.time = 0;
    ave.label = CM.label;
    ave.dimord = 'chan_time';
    ft_topoplotER(cfg,ave);
end









% Assessing statistical significance

% extracting CKC and plotting

subjects = subjects_kids; group = 'kids';
subjects = subjects_adults; group = 'adults';
subjects = cat(2,subjects_kids,subjects_adults); group = 'both';


refs = {'vCOP' 'abc_vCOP' 'rCoP'};
p = [];
for n_sub = 1:length(subjects)
    for n_cond = 1:4
        mat_outfile = fullfile(matfold,[subjects{n_sub} '_' conds{n_cond} '_coherence.mat']);
        load(mat_outfile,'surrog','CMall')
        surrog_computed(n_sub,n_cond) = isempty(surrog);

        for n_ref = 1:3

            if 1 % same electrode selection as in Thomas'
                if not(exist('SM1'))
                    to_include = {'Cz' 'C1' 'C2' 'FCz' 'FC1' 'FC2' 'CP1' 'CP2'};
                    SM1 = [];
                    for n = 1:length(to_include)
                        SM1 = [SM1 find(strcmp(CMall.label,to_include{n}))];
                    end
                end
                keepE = SM1;
            end

            keepF = find(CMall.f <= 6 & CMall.f >= 1);

            Coh_genuine = max(max(CMall.coh(keepE,keepF,n_ref)));
            Coh_surrog = squeeze(max(max(surrog.coh(keepE,keepF,n_ref,:))));
            p(n_sub,n_cond,n_ref) = mean(Coh_genuine <= Coh_surrog);
        end
    end
end

squeeze(sum(p <= 0.05))  % number of significant subjects

squeeze(sum(p <= 0.05))/size(p,1)*100  % percentage of significant subjects





% extracting CKC and plotting

subjects = subjects_kids; group = 'kids';
subjects = subjects_adults; group = 'adults';
subjects = cat(2,subjects_kids,subjects_adults); group = 'both';

surflap = 0;
make_plots = 0;
ckcMean = [];
nave = [];
refs = {'vCOP', 'abc_vCOP', 'rCoP'};

for n_ref = 1:3
    group_CKC_trace = [];
    group_CKC_topo = [];
    for n_sub = 1:length(subjects)
        if make_plots
            dim = get(0, 'ScreenSize');
            fig = figure('Menubar', 'none', 'Visible', 'Off');
            set(fig, 'Units', 'pixel', 'Position', dim)
        end
        for n_cond = 1:4
            mat_outfile = fullfile(matfold, [subjects{n_sub} '_' conds{n_cond} '_coherence.mat']);
            load(mat_outfile, 'CMall')
            CMall.Fyy(1,1:2,n_ref) = 0;

            nave(n_sub, n_cond) = CMall.nave;

            if ~exist('keepE', 'var')
                outer_rim = {'M1', 'M2', 'Oz', 'O2', 'PO8', 'P8', 'TP8', 'T8', 'FT8', 'F8', 'AF8', 'Fp2', 'Fpz', 'Fp1', 'AF7', 'F7', 'FT7', 'T7', 'TP7', 'P7', 'PO7', 'O1'};
                to_rm = [];
                for n = 1:length(outer_rim)
                    to_rm = [to_rm find(strcmp(CMall.label, outer_rim{n}))];
                end
                keepE = 1:61;
                keepE(to_rm) = [];
            end

            to_include = {'Cz', 'C1', 'C2', 'FCz', 'FC1', 'FC2', 'CP1', 'CP2'};
            SM1 = [];
            for n = 1:length(to_include)
                SM1 = [SM1 find(strcmp(CMall.label, to_include{n}))];
            end
            keepE = SM1;

            coh_smooth = CMall.coh(keepE, :, n_ref);
            keepF = find(CMall.f >= 1 & CMall.f <= 6);
            M = max(coh_smooth(:, keepF), [], 2);
            [M, thisElecIdx] = max(M);
            thisElec = keepE(thisElecIdx);

            ckcMean(n_sub, n_cond, n_ref) = M;
            group_CKC_trace(:, n_sub, n_cond) = CMall.coh(thisElec, :, n_ref);

            ave = [];
            ave.avg = mean(CMall.coh(1:61, keepF, n_ref), 2);
            group_CKC_topo(:, n_sub, n_cond) = ave.avg;

            if make_plots
                subplot(4,3,(1:2)+(n_cond-1)*3)
                plot(CMall.f, CMall.coh(keepE, :, n_ref)', 'k')
                hold on
                plot(CMall.f, CMall.coh(thisElec, :, n_ref), 'r', 'linewidth', 2)
                ylim([0 0.16])  % range fisso asse y
                if n_cond == 4
                    xlabel('Frequency [Hz]')
                end
                ylabel('Coherence')
                ttt = conds{n_cond}; ttt(ttt == '_') = ' ';
                title(ttt)

                cfg = [];
                cfg.zlim = [0 0.008]; % range fisso topoplot
                cfg.comment = 'no';
                cfg.layout = 'EEG1010.lay';
                ave.fsample = 1000;
                ave.time = 0;
                ave.label = CMall.label(1:61);
                ave.dimord = 'chan_time';
                subplot(4,3,3+(n_cond-1)*3)
                ft_topoplotER(cfg, ave);
                colorbar
            end
        end
        if make_plots
            picsname = fullfile(picsfold, ['CKC_' subjects{n_sub} '_' refs{n_ref} '.jpeg']);
            print(picsname, '-djpeg');
            close(fig)
        end
    end

    dim = get(0, 'ScreenSize');
    fig = figure('Menubar', 'none', 'Visible', 'On');
    set(fig, 'Units', 'pixel', 'Position', dim)
    for n_cond = 1:4
        % subplot(4,3,(1:2)+(n_cond-1)*3)
        axes('Units','Centimeters','Position',[2 22.5-n_cond*5 5 4])
        plot(CMall.f, group_CKC_trace(:, :, n_cond), 'k')
        ylim([0 0.16])  % range fisso asse y
        if n_cond == 4
            xlabel('Frequency [Hz]')
        end
        ylabel('Coherence')
        ttt = conds{n_cond}; ttt(ttt == '_') = ' ';
        title(ttt)

        ave = [];
        ave.avg = mean(group_CKC_topo(:, :, n_cond), 2);
        cfg = [];
        cfg.zlim = [0 0.008]; % range fisso topoplot
        cfg.comment = 'no';
        cfg.layout = 'EEG1010.lay';
        ave.fsample = 1000;
        ave.time = 0;
        ave.label = CMall.label(1:61);
        ave.dimord = 'chan_time';
        % subplot(4,3,3+(n_cond-1)*3)
        axes('Units','Centimeters','Position',[8 22.5-n_cond*5 4 4])
        ft_topoplotER(cfg, ave);
        colorbar
    end
    picsname = fullfile(picsfold, ['CKC_' group '_' refs{n_ref} '.jpeg']);
    print(picsname, '-djpeg');
    picsname = fullfile(picsfold, ['CKC_' group '_' refs{n_ref} '.svg']);
    print(picsname, '-dsvg');
    close(fig)
end


n_ref = 1;
CM_anovan(ckcMean(:,:,n_ref),{'subjects' 'cond' 'group'},[length(subjects_kids) length(subjects_adults)]);
 

CM_anovan(nave,{'subjects' 'cond' 'group'},[length(subjects_kids) length(subjects_adults)]);


%load nave nave

nave1 = nave(1:length(subjects_kids),:);
nave2 = nave(length(subjects_kids)+1:end,:);
nave2 = nave(length(subjects_kids)+1:2*length(subjects_kids),:);

nave_low = min(mean(nave1))
nave1 = nave1 - (mean(nave1)-nave_low);

nave2old = nave2;
for n_cond = 1:4
    [val1,order1] = sort(nave1(:,n_cond));
    [val2,order2] = sort(nave2(:,n_cond));
    nave2(order2,n_cond) = round(nave1(order1,n_cond));

    while any(nave2(:,n_cond) > nave2old(:,n_cond))
        n_to_perm = find(nave2(:,n_cond) > nave2old(:,n_cond));
        n_to_perm_with = 1:length(nave2);
        n_to_perm_with(n_to_perm) = [];
        n_to_perm_with = n_to_perm_with(randperm(length(n_to_perm_with),length(n_to_perm)));

        vals = nave2(n_to_perm,n_cond)-nave2old(n_to_perm,n_cond)
        nave2(n_to_perm,n_cond) = nave2(n_to_perm,n_cond)-vals;
        nave2(n_to_perm_with,n_cond) = nave2(n_to_perm_with,n_cond)+vals;
    end
end

nave = cat(1,round(nave1),nave2,nave(2*length(subjects_kids)+1:end,:))


save nave_corrected nave subjects




% Media e deviazione standard globali (across all participants and conditions)
m_epochs  = mean(nave(:));
sd_epochs = std(nave(:));

fprintf('Across participants and conditions, %.0f ± %.0f epochs were retained for CKC estimation\n', m_epochs, sd_epochs);















%%
% for visualization
n_sub = 1;
n_cond = 4;
mat_outfile = fullfile(matfold,subjects{n_sub},[subjects{n_sub} '_' conds{n_cond} '_coherence.mat']);
load(mat_outfile)

Neeg = 61;
n_ref = 2; % 1 biceps non-rectif; 2 biceps rectif; 3 FDI non-rectif; 4 FDI rectif
F = 3.5;
Find = round(F(1)/CMall.Fs*CMall.quantum+1):round(F(end)/CMall.Fs*CMall.quantum+1);
subplot(1,2,1)
plot(CMall.f,CMall.coh(1:Neeg,:,n_ref)')
xlabel('Frequency [Hz]'); ylabel('Coherence')
if 1
    % colmap = colormap('jet');

    cfg = [];
    cfg.xparam = 'time';
    % cfg.zlim = [0 0.16]
    cfg.layout = 'EEG1010.lay';
    ave = [];
    ave.avg = mean(CMall.coh(1:Neeg,Find,n_ref),2);
    cfg.zlim = [0 max(ave.avg)];
    cfg.comment = 'no';
    ave.fsample = 1000;
    ave.time = 0;
    ave.label = CMall.label(1:Neeg);
    ave.dimord = 'chan_time';
    subplot(1,2,2)
    ft_topoplotER(cfg,ave);
    title(['Coherence at ' num2str(F) ' Hz'])
end

F = [3 5];
Find = round(F(1)/CMall.Fs*CMall.quantum+1):round(F(end)/CMall.Fs*CMall.quantum+1);
[M,n_elec] = max(mean(CMall.coh(1:Neeg,Find,n_ref),2));

subplot(2,2,1)
plot(CMall.f,(CMall.coh(n_elec,:,n_ref)))
xlabel('Frequency [Hz]'); ylabel('Coherence')

subplot(2,2,3)
plot(CMall.f,phase(CMall.Fxy(n_elec,:,n_ref)),'o')
xlabel('Frequency [Hz]'); ylabel('Phase')

subplot(2,2,2)
plot((-CMall.quantum/2+1:CMall.quantum/2)/Fs*1000,CMall.fxy(n_elec,:,n_ref))
xlabel('Time [ms]'); title('cross-correlogram')




plot (CMall.fxy(n_elec,:,n_ref)


plot(CMderivall.f,CMderivall.coh(1:Nderiv,:,n_ref)')

to_plot_deriv = CMderivall.coh(1:Nderiv,find(CMall.f==F),n_ref);
to_plot = zeros(Neeg,1);
for n_elec = 1:Neeg
    in_deriv = find(any(all_derivs' == n_elec));
    to_plot(n_elec) = mean(to_plot_deriv(in_deriv));
end


if 0
    cfg = [];
    cfg.xparam = 'time';
    % cfg.zlim = [0 0.16]
    cfg.layout = 'EEG1010.lay';
    ave = [];
    ave.avg = to_plot;
    cfg.zlim = [0 max(ave.avg)];
    cfg.comment = 'no';
    ave.fsample = 1000;
    ave.time = 0;
    ave.label = CMall.label;
    ave.dimord = 'chan_time';
    ft_topoplotER(cfg,ave);
end
