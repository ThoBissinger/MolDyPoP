%% initialization
clear
cd /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/EquilibrationPlots/Scaling
run('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/initialization_script.m')
pathbase=pwd;
model="mxy";
sqrtN = 64;
% sqrtN_vals = [64];
step_str="16";
% step_str=["2" "4" "8"];
step=str2double(step_str);
T_str=[".14", ".185", ".205"];
T_vals=str2double(T_str);
N_T=numel(T_vals);

filebase = "/data/scc/thobi/220829_EquilibrationCheck/anneal/ann_step_";
filemid = "/aligned/mxy_3.00/sqrtN_";
fileend = "/T_.01/samp_eq";
filename = filebase + step_str + filemid + sqrtN + fileend;
if step == 2
    filename = "/data/scc/thobi/220829_EquilibrationCheck/anneal" + filemid + sqrtN + fileend;
end
S_long=load(filename);

for i_T = 1:N_T
    filename="/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_64/T_" + T_str(i_T) + "/eq_samp_eq";
    S{i_T}=load(filename);
end

%% m vs temperature dependence
c_map = lines(N_T);
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);

T=S_long.temperature(2:end);
m=S_long.absM(2:end);
plot(T,m,'-','Color','k',...
    'DisplayName',sprintf('$\\tau = %s$',step_str),...
    "LineWidth",1.3);
hold on;
for i_T = 1:N_T
    T_target = T_vals(i_T);
    T = S{i_T}.temperature(1:end);
    m = S{i_T}.absM(1:end);
    dispname = sprintf('$T = %.3f$',T_target);
%         dispname = sprintf('$\\tau_{\\textrm{anneal}} = %.2g$',step_vals(i_step));
    plot(T,m,'DisplayName',dispname,...
        'Color',c_map(i_T,:),...
        "LineWidth",1.3);
    hold on;
end
xlim([0 .4]);
hlegend=legend('Location','SouthWest','Interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_annotation);
%     legpos=get(hlegend,"Position");
%     legpos=legpos+[0 .1 0 0];
%     set(hlegend,"Position",legpos)
hXLabel = xlabel('$\langle T\rangle$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$\langle m \rangle$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

%     text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');

figname=sprintf('%s/%s_sqrtN_%d_m_vs_temp',pathbase,model,sqrtN);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
