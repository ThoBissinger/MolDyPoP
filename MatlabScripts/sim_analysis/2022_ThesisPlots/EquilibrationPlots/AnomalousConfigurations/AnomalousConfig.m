clear
cd /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/EquilibrationPlots/AnomalousConfigurations
run('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/initialization_script.m')
pathbase=pwd;
model="mxy";
sqrtN = 128;
S=load('/data/scc/thobi/220829_EquilibrationCheck/anneal/aligned/mxy_3.00/sqrtN_128/T_.01/samp_eq_collect',"absM_collect","averaging_times");

tt = S.averaging_times;
mm = S.absM_collect;
img_res=600;

t_vals = [2e3 4e3 6e3 8e3 1e4];
N_t = numel(t_vals);
%% Energy vs t
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);
c_map = lines(N_t);
for i_t = 1:N_t
    t = t_vals(i_t);
    [minval,index] = min(abs(tt - t));
    dispname = sprintf('$t = %d$',t);
    plot(1:100,mm(:,index),'--o','DisplayName',dispname,...
        'Color',c_map(i_t,:));
    hold on;
end
hlegend=legend('Location','SouthEast','Interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_annotation,"NumColumns",1);
hXLabel = xlabel('run number','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$m(t)$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

% text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');

figname=sprintf('%s/%s_%d_m_runs',pathbase,model,sqrtN);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',img_res);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');


%% Snapshot of anomalous configuration
img_size=[0 0 2*columnwidth_cm 2*columnwidth_cm];
figure
set(gcf,'units','centimeters','OuterPosition',img_size);
single_snap_fig("/data/scc/thobi/220829_EquilibrationCheck/anneal/aligned/mxy_3.00/sqrtN_128/T_.01/run_33/output/snapshot_Dynamics_final.out",'mxy','spins',"s");
xlim([0 74]);
ylim([0 74]);
axis off
% axis tight
% hXLabel = xlabel('x','interpreter','latex', ...
%     'FontName', 'cmr12','FontSize', fontsize_ax_labels);
% hYLabel = ylabel('$y$','interpreter','latex', ...
%     'FontName', 'cmr12','FontSize', fontsize_ax_labels);
figname=sprintf('%s/%s_%d_run_33',pathbase,model,sqrtN);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',img_res);

figure
set(gcf,'units','centimeters','OuterPosition',img_size);
single_snap_fig("/data/scc/thobi/220829_EquilibrationCheck/anneal/aligned/mxy_3.00/sqrtN_128/T_.01/run_11/output/snapshot_Dynamics_final.out",'mxy','spins',"s");
xlim([0 74]);
ylim([0 74]);
axis off
% axis tight
figname=sprintf('%s/%s_%d_run_11',pathbase,model,sqrtN);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',img_res);



figure
set(gcf,'units','centimeters','OuterPosition',img_size);
single_snap_fig("/data/scc/thobi/220829_EquilibrationCheck/anneal/aligned/mxy_3.00/sqrtN_128/T_.01/run_30/output/snapshot_Dynamics_final.out",'mxy','spins',"s");
xlim([0 74]);
ylim([0 74]);
axis off
% axis tight
figname=sprintf('%s/%s_%d_run_30',pathbase,model,sqrtN);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',img_res);


figure
set(gcf,'units','centimeters','OuterPosition',img_size);
single_snap_fig("/data/scc/thobi/220829_EquilibrationCheck/anneal/aligned/mxy_3.00/sqrtN_128/T_.01/run_14/output/snapshot_Dynamics_final.out",'mxy','spins',"s");
xlim([0 74]);
ylim([0 74]);
axis off
% axis tight
figname=sprintf('%s/%s_%d_run_14',pathbase,model,sqrtN);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',img_res);
