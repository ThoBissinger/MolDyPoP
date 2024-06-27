clear
cd /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/EquilibrationPlots/TimeDependence
run('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/initialization_script.m')
pathbase=pwd;
model="mxy";
sqrtN_vals = [16, 32, 64, 128];
N_N=numel(sqrtN_vals);

srcdir="/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00";
% T_str=[".01" ".03" ".05" ".07" ".09" ".11" ".13" ".14" ".15" ".155" ".16" ".165" ".17" ".175" ".18" ".185" ".19" ".195" ".20" ".205" ".21" ".22" ".23" ".24" ".25" ".27" ".29" ".31" ".33" ".35" ".37" ".40" ".43" ".46" ".49" ".52"];
T_str=[".01" ".03" ".05" ".07" ".09" ".11" ".13" ".14" ".15" ".16" ".17" ".18" ".185" ".19" ".195" ".20" ".205" ".21" ".22" ".23" ".24" ".25" ".27" ".29" ".31" ".33" ".35" ".37"];
T_vals=str2double(T_str);
N_T = numel(T_str);
for i_N = 1:N_N
    sqrtN=sqrtN_vals(i_N);
    S{i_N}=load("/data/scc/thobi/220829_EquilibrationCheck/anneal/aligned/mxy_3.00/sqrtN_" ...
        + sqrtN + "/T_.01/samp_eq");
    for i_T = 1:N_T
        S_aux=load(srcdir + "/sqrtN_" + sqrtN + "/T_" + T_str(i_T) + "/samp_Dynamics","absM_av");
        absM_vals_aux(i_T) = S_aux.absM_av;
    end
    absM_cell{i_N} = absM_vals_aux;
end

%% Energy vs t
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);
c_map = linspecer(N_N);
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    t = S{i_N}.averaging_times;
    E = S{i_N}.H;
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    plot(t,E,'DisplayName',dispname,...
        'Color',c_map(i_N,:));
    hold on;
end
hlegend=legend('Location','SouthEast','Interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_annotation);
legpos=get(hlegend,"Position");
legpos=legpos+[0 .1 0 0];
set(hlegend,"Position",legpos)
hXLabel = xlabel('$t$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$\langle H \rangle/ N$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');

ax_inset = axes('Position',[.5 .5 .38 .38],'Units','Normalized');
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    t = S{i_N}.averaging_times;
    E = S{i_N}.H;
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    plot(t,E,'-','DisplayName',dispname,...
        'Color',c_map(i_N,:), 'MarkerSize',4);
    hold on;
end
xlim([0 30]);

figname=sprintf('%s/%s_E_vs_t',pathbase,model);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');



%% Tempearture vs t
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);
c_map = linspecer(N_N);
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    t = S{i_N}.averaging_times;
    T = S{i_N}.temperature;
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    plot(t,T,'DisplayName',dispname,...
        'Color',c_map(i_N,:));
    hold on;
end
% hlegend=legend('Location','SouthEast','Interpreter','latex', ...
%     'FontName', 'cmr12','FontSize', fontsize_annotation);
% legpos=get(hlegend,"Position");
% legpos=legpos+[0 .1 0 0];
% set(hlegend,"Position",legpos)
hXLabel = xlabel('$t$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$\langle T\rangle$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

text(0.025,0.95,"b)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');

ax_inset = axes('Position',[.5 .5 .38 .38],'Units','Normalized');
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    t = S{i_N}.averaging_times;
    T = S{i_N}.temperature;
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    plot(t,T,'DisplayName',dispname,...
        'Color',c_map(i_N,:));
    hold on;
end
xlim([0 30]);
ylim([.46 .52]);

figname=sprintf('%s/%s_T_vs_t',pathbase,model);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');




%% M vs t
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);
c_map = linspecer(N_N);
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    t = S{i_N}.averaging_times;
    M = S{i_N}.absM;
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    plot(t,M,'DisplayName',dispname,...
        'Color',c_map(i_N,:));
    hold on;
end
% hlegend=legend('Location','SouthEast','Interpreter','latex', ...
%     'FontName', 'cmr12','FontSize', fontsize_annotation);
% legpos=get(hlegend,"Position");
% legpos=legpos+[0 .1 0 0];
% set(hlegend,"Position",legpos)
hXLabel = xlabel('$t$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$\langle m\rangle$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');

ax_inset = axes('Position',[.5 .18 .38 .5],'Units','Normalized');
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    t = S{i_N}.averaging_times;
    M = S{i_N}.absM;
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    plot(t,M,'DisplayName',dispname,...
        'Color',c_map(i_N,:));
    hold on;
end
xlim([0 80]);

figname=sprintf('%s/%s_m_vs_time',pathbase,model);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');



%% M vs T
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);
c_map = linspecer(N_N);


for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    T = S{i_N}.temperature(2:end);
    M = S{i_N}.absM(2:end);
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    plot(T,M,'-','DisplayName',dispname,...
        'Color',c_map(i_N,:),'LineWidth',1.5);
    hold on;
end
for i_N = 1:N_N
    plot(T_vals,absM_cell{i_N},'--','DisplayName',dispname,...
        'Color','k','LineWidth',1);
    hold on;
end
xlim([0 .38])
% hlegend=legend('Location','SouthEast','Interpreter','latex', ...
%     'FontName', 'cmr12','FontSize', fontsize_annotation);
% legpos=get(hlegend,"Position");
% legpos=legpos+[0 .1 0 0];
% set(hlegend,"Position",legpos)
hXLabel = xlabel('$\langle T \rangle$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$\langle m\rangle$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

text(0.025,0.05,"b)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');



figname=sprintf('%s/%s_m_vs_Temp',pathbase,model);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
