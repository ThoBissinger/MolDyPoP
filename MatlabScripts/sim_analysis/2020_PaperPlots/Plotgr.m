clear all
close all
% addpath /home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;

mxydata=load('mxy/rho_3.00_grselect.mat'); 

dirs=mxydata.dirs;
rbin=mxydata.currbin;
gr=mxydata.gr;

indices = find(rbin < 3);
rbin=rbin(indices);

T_str=[".01" ".09" ".15" ".19" ".25" ".31"];

N_files = numel(dirs);
for i = 1 : N_files
    curgr{i}= gr{i,1}(indices);
    nruns=1;
    for runnr=2:125
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
for i = 1 : N_files
    index=i;
    hgrplot{i} = plot(rbin,curgr{index});
    dispname{i}=sprintf('$T = %s$', T_str(index));
    set(hgrplot{i}, ...
        'LineStyle', '-', 'LineWidth', 2, ...
        'DisplayName', dispname{i}, ...
        'Color', c_map(index,:))
end
ylim([0 1.2*max(curgr{1})]);
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
for i = 1 : N_files-1
    index=i+1;
    hgrplot{i} = plot(rbin,curgr{index});
    dispname{i}=sprintf('$T = %s$', T_str(index));
    set(hgrplot{i}, ...
        'LineStyle', '-', 'LineWidth', 2, ...
        'DisplayName', dispname{i}, ...
        'Color', c_map(index,:))
end
xlim([.35 1.8]);
ylim([.8 1.4]);
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


figname='plots/mxy_gr';
fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end
