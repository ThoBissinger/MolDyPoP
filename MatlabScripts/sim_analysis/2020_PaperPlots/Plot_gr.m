clear all
close all

initialization_script;

basedir=sprintf('%s/plots/gr',fig_base);

saveswitch=1;

model = "mxy";
% model = "NoSpin";
if ( model == "mxy")
    data=load('mxy/rho_3.00_grselect.mat'); 
    T_str=[".01" ".03" ".09" ".15" ".19" ".31"]; % ".25" 
    runmax=500;
    figname='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/mxy_gr';
elseif ( model == "NoSpin")
    data=load('NoSpin/NoSpin_grselect.mat'); 
    T_select=[1,5,8,10,12,13];
    T_vals=[.01, ];
    N_files=numel(T_select);
    runmax=50;
    figname='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/NoSpin_gr';
end


% sqrtN_vals = [16, 32, 64, 128, 256];
sqrtN_vals = [16,32,64,128];
sqrtN_dirs = sprintfc('/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_%d',sqrtN_vals);
N_N = numel(sqrtN_vals);
N_T = numel(T_str);

c_map=linspecer(N_T);
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    figure
    hold on;
    dispname = sprintf('$T = %s$',T_str{1});
    plot(nan,nan,'-',...
        'LineStyle', '-', 'LineWidth', 1.2, ...
        'DisplayName', dispname, ...
        'Color', c_map(1,:));
    for i_T = 2:N_T
        curdir = sprintf('%s/T_%s/',sqrtN_dirs{i_N},T_str{i_T});
        curfile = sprintf('%s/samp_Dynamics_gr',curdir);
        dispname = sprintf('$T = %s$',T_str{i_T});
        load(curfile,"gr","rbin");
        plot(rbin,gr,'-',...
            'LineStyle', '-', 'LineWidth', 1.2, ...
            'DisplayName', dispname, ...
            'Color', c_map(i_T,:));
        
    end
    xlim([.35 1.8]);
    ylim([0 2.1]);
    hLegend=legend('Location','southeast','Interpreter','latex',...
        'NumColumns',2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Adjust axes properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')
    
    hXLabel = xlabel('$r$','interpreter','latex');
    hYLabel = ylabel('$g(r)$','interpreter','latex');
    h_axis = gca;
    
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    
    set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
    
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inset
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax_inset = axes('Position',[.59 .61 .3 .3]);
    
    set(ax_inset, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'LineWidth', .5, 'Xscale', 'lin')
    maxgr=0;
    for i_T = 1 : N_T
        curdir = sprintf('%s/T_%s/',sqrtN_dirs{i_N},T_str{i_T});
        curfile = sprintf('%s/samp_Dynamics_gr',curdir);
        dispname = sprintf('$T = %s$',T_str{i_T});
        load(curfile,"gr","rbin");
        plot(rbin,gr,'-',...
            'LineStyle', '-', 'LineWidth', 1.2, ...
            'DisplayName', dispname, ...
            'Color', c_map(i_T,:));
        hold on;
        maxgr = max(maxgr,max(gr));
    end
    ylim([0 1.2*maxgr]);
    xlim([0 4]);
%     [ymax,i_max] = max(curgr{1});
%     text(rbin(i_max),ymax,'$T = .01$',...
%         'VerticalAlignment','bottom','HorizontalAlignment','left',...
%         'interpreter','latex',...
%         'Color','red','FontSize', 10);
    
    % set(hInsetPlot, ...
    %     'LineStyle', '--', 'LineWidth', 1.5, ...
    %     'Marker', '+', 'MarkerSize', 6, ...
    %     'Color',c_map(1,:))
    
    %     hXLabel = xlabel('$T$','interpreter','latex');
    %     hYLabel = ylabel('$D_r$','interpreter','latex');
    NW = [min(xlim) max(ylim)]+[1.6*diff(xlim) -diff(ylim)]*0.02;
    SE = [max(xlim) min(ylim)]+[-diff(xlim) .5*diff(ylim)]*0.02;
    text(SE(1),SE(2),'$r$',...
        'VerticalAlignment','bottom','HorizontalAlignment','right',...
        'interpreter','latex',...
        'Color','black','FontSize', fontsize_ax_labels);
    text(NW(1),NW(2),'$g(r)$',...
        'VerticalAlignment','top','HorizontalAlignment','left',...
        'interpreter','latex',...
        'Color','black','FontSize', fontsize_ax_labels);
    set(gca, 'FontName', 'cmr12')
    %     set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set([gca], 'FontSize', fontsize_axis)
    

    figname=sprintf('%s/mxy_gr_sqrtN_%d',basedir,sqrtN);
    
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end

return
% runmax=125; %125

dirs=data.dirs;
fullrbin=data.rbin;
rbin=fullrbin{end,end};
gr=data.gr;
T_str=data.T_str;
T_vals=data.T_vals;

indices = find(rbin < 3);
rbin=rbin(indices);

for i = 1 : numel(dirs)
    
    curgr{i}= gr{i,1}(indices);
    nruns=1;
    for runnr=2:runmax
        if (~ isempty(gr{i,runnr}))
            curgr{i} = curgr{i} + gr{i,runnr}(indices);
            nruns=nruns+1;
        end
    end
    curgr{i} = curgr{i}/nruns;
end
% mxyfit=load('mxy/rho_3.00_CritExpFit.mat');
% xydata=load('xy/lf0_eq.mat'); xyfit=load('xy/lf0_CritExpFit.mat');


% files = ["/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_256/T_.01/run_1/output/snapshot_eq_final.out" ...
%     "/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_256/T_.15/run_1/output/snapshot_eq_final.out" ...
%     "/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_256/T_.19/run_1/output/snapshot_eq_final.out" ...
%     "/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_256/T_.25/run_1/output/snapshot_eq_final.out"];


% files = ["/data/scc/thobi/201207_equilibration/xy/anneal/sqrtN_64/T_.10/run_1/output/snapshot_eq_final.out" ...
%     "/data/scc/thobi/201207_equilibration/xy/anneal/sqrtN_64/T_1.30/run_1/output/snapshot_eq_final.out" ...
%     "/data/scc/thobi/201207_equilibration/xy/anneal/sqrtN_64/T_1.50/run_1/output/snapshot_eq_final.out" ...
%     "/data/scc/thobi/201207_equilibration/xy/anneal/sqrtN_64/T_1.60/run_1/output/snapshot_eq_final.out"];


figure(1)
hold on;
c_map=linspecer(N_files);
for i = 1 : N_files-1
    index=T_select(i+1);
    hgrplot{i} = plot(rbin,curgr{index});
    dispname{i}=sprintf('$T = %s$', T_str(index))
    set(hgrplot{i}, ...
        'LineStyle', '-', 'LineWidth', 1.2, ...
        'DisplayName', dispname{i}, ...
        'Color', c_map(i+1,:))
end
xlim([.35 1.8]);
ylim([.8 1.5]);
legend('Location','northwest','Interpreter','latex');
%% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')

hXLabel = xlabel('$r$','interpreter','latex');
hYLabel = ylabel('$g(r)$','interpreter','latex');
ax_full = gca;


    
%% Inset
ax_inset = axes('Position',[.53 .53 .35 .35]);
hold on;
set(ax_inset, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', .5, 'Xscale', 'lin')
for i = 1 : N_files
    index=T_select(i);
    hgrplot_inset{i} = plot(rbin,curgr{index});
    dispname{i}=sprintf('$T = %s$', T_str(index));
    set(hgrplot_inset{i}, ...
        'LineStyle', '-', 'LineWidth', 2, ...
        'DisplayName', dispname{i}, ...
        'Color', c_map(i,:))
end
ylim([0 1.2*max(curgr{1})]);
[ymax,i_max] = max(curgr{1});
text(rbin(i_max),ymax,'$T = .01$',...
    'VerticalAlignment','bottom','HorizontalAlignment','left',...
    'interpreter','latex',...
    'Color','red','FontSize', 10);

% set(hInsetPlot, ...
%     'LineStyle', '--', 'LineWidth', 1.5, ...
%     'Marker', '+', 'MarkerSize', 6, ...
%     'Color',c_map(1,:))

%     hXLabel = xlabel('$T$','interpreter','latex');
%     hYLabel = ylabel('$D_r$','interpreter','latex');
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),'$r$',...
    'VerticalAlignment','bottom','HorizontalAlignment','right',...
    'interpreter','latex',...
    'Color','black','FontSize', 10);
text(NW(1),NW(2),'$g(r)$',...
    'VerticalAlignment','top','HorizontalAlignment','left',...
    'interpreter','latex',...
    'Color','black','FontSize', 10);
set(gca, 'FontName', 'cmr12')
%     set([hXLabel, hYLabel], 'FontName', 'cmr12')
set([gca], 'FontSize', 8)
%     set([hXLabel, hYLabel], 'FontSize', 12)


% figname='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/mxy_gr';
fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end
