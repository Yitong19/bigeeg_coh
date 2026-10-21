[max_v_ref,loc_v_ref] = findpeaks(squeeze(CM.Fyy(:,:,1)));  
[~, order_v_ref] = sort(max_v_ref, 'descend');
[max_y_ref,loc_y_ref] = findpeaks(squeeze(CM.Fyy(:,:,2))); 
[~, order_y_ref] = sort(max_y_ref, 'descend');


[max_v,loc_v] = findpeaks(squeeze(mean(CM.coh(1:64,:,1))));  
[~, order_v] = sort(max_v, 'descend');

[max_y,loc_y] = findpeaks(squeeze(mean(CM.coh(1:64,:,2)))); 
[~, order_y] = sort(max_y, 'descend');

figure('Units','normalized', 'Position',[0.05 0.12 0.90 0.72]);

t = tiledlayout(2,4, ...
    'TileSpacing','compact', ...
    'Padding','compact');

nexttile(1);
hData_v = plot(CM.f, squeeze(CM.Fyy(:,:,1))');
set(hData_v, 'HandleVisibility', 'off');  
hold on;

f1_v = CM.f(loc_v_ref(order_v_ref(1)));
f2_v = CM.f(loc_v_ref(order_v_ref(2)));

h1_v = xline(f1_v, '--r', 'LineWidth', 1.5);
h2_v = xline(f2_v, '--b', 'LineWidth', 1.5);

legend([h1_v h2_v], ...
    {sprintf('%.2f Hz', f1_v), ...
     sprintf('%.2f Hz', f2_v)}, ...
    'Location', 'best');

hold off;

nexttile(5);
hData_y = plot(CM.f, squeeze(CM.Fyy(:,:,2))');
set(hData_y, 'HandleVisibility', 'off');  
hold on;

f1_y = CM.f(loc_y_ref(order_y_ref(1)));

h1_y = xline(f1_y, '--r', 'LineWidth', 1.5);

legend(h1_y,sprintf('%.2f Hz', f1_y),'Location', 'best');


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


nexttile(6);
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

ave_v0 = [];
ave_v0.avg     = squeeze(CM.coh(1:64,loc_v(order_v(1)),1));
ave_v0.fsample = 1000;
ave_v0.time    = 0;
ave_v0.label   = data_cm.label(1:64);
ave_v0.dimord  = 'chan_time';

nexttile(3);
ft_topoplotER(cfg, ave_v0);
title(sprintf('%.2f Hz', CM.f(loc_v(order_v(1)))));

ave_v1 = [];
ave_v1.avg     = squeeze(CM.coh(1:64,loc_v(order_v(2)),1));
ave_v1.fsample = 1000;
ave_v1.time    = 0;
ave_v1.label   = data_cm.label(1:64);
ave_v1.dimord  = 'chan_time';

nexttile(4);
ft_topoplotER(cfg, ave_v1);
title(sprintf('%.2f Hz', CM.f(loc_v(order_v(2)))));

ave_y0 = [];
ave_y0.avg     = squeeze(CM.coh(1:64,loc_y(order_y(1)),2));
ave_y0.fsample = 1000;
ave_y0.time    = 0;
ave_y0.label   = data_cm.label(1:64);
ave_y0.dimord  = 'chan_time';

nexttile(7);
ft_topoplotER(cfg, ave_y0);
title(sprintf('%.2f Hz', CM.f(loc_y(order_y(1)))));

ave_y1 = [];
ave_y1.avg     = squeeze(CM.coh(1:64,loc_y(order_y(2)),2));
ave_y1.fsample = 1000;
ave_y1.time    = 0;
ave_y1.label   = data_cm.label(1:64);
ave_y1.dimord  = 'chan_time';

nexttile(8);
ft_topoplotER(cfg, ave_y1);
title(sprintf('%.2f Hz', CM.f(loc_y(order_y(2)))));


xy = cat(1,v_all_cat,y_all_cat,ts_all_cat);
           
            CM = [];
            CM.ds = 1;
            CM.p = p;
            CM.Fmax = 10;
            CM.Fs = Fs;
            CM.to_compute.Coh = 1;
            CM.to_compute.rPDC = 1;
            CM.show_plots = 0;
            CM.stat.alpha = 0.05;
            % CM.smooth = 3;
            CM = CM_rPDC_stab(CM,xy,bad);
            rPDC(1,:,n_s,n_sub,n_cond) = squeeze(CM.rPDC(1,2,:));
            rPDC(2,:,n_s,n_sub,n_cond) = squeeze(CM.rPDC(2,1,:));
            Coh(:,n_s,n_sub,n_cond) = squeeze(CM.Coh(1,2,:));
            phi(1,:,n_s,n_sub,n_cond) = unwrap(angle(squeeze(CM.FA(1,2,:))));
            phi(2,:,n_s,n_sub,n_cond) = unwrap(angle(squeeze(CM.FA(2,1,:))));
            
            % time-reversed rPDC
            CM = CM_rPDC_stab(CM,xy([n_s 5],end:-1:1),xy(6,end:-1:1));
            rPDC_rev(1,:,n_s,n_sub,n_cond) = squeeze(CM.rPDC(1,2,:));
            rPDC_rev(2,:,n_s,n_sub,n_cond) = squeeze(CM.rPDC(2,1,:));