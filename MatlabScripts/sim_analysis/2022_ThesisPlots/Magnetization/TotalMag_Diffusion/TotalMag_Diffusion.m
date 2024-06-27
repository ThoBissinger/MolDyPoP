%% initialization
clear
pathbase="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/Magnetization/TotalMag_Diffusion";
cd(pathbase)
run('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/initialization_script.m')
% pathbase=pwd;
model="mxy";


if model == "mxy"
    sqrtN_vals = [16,32,64,128,256];
    sqrtN = 16;
    T_str = [".01" ".03" ".05" ".07" ".09" ".11" ".13" ".14" ".15" ".155" ".16" ".165" ".17" ".175" ".18" ".185" ".19" ".195" ".20" ".205" ".21" ".22" ".23" ".24" ".25" ".27" ".29" ".31" ".33" ".35" ".37" ".40" ".43" ".46" ".49" ".52" ];
    T_str = [".07" ".11" ".14" ".17" ".19" ".20"];
    filebase = "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_";
    fileend="/samp_Dynamics_M_TimeEvolution.mat";
    filebase_shorttime = "/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00/sqrtN_";
    filename="samp_Dynamics_M_TimeEvolution";
elseif model == "xy"
    sqrtN_vals = [16,32,64,128];
%     T_str = [".10" ".20" ".30" ".40" ".50" ".60" ".70" ".80" ".85" ".91" ".95" "1.00"];
    T_str = [".10" ".20" ".30" ".40" ".50" ".60" ".70" ".80" ".85" ".87" ".89" ".90" ".91" ".93" ".95" ".97" "1.00" "1.03" "1.06" "1.09" "1.10" "1.12" "1.15" "1.18" "1.20" "1.21" "1.24" "1.30" "1.40" "1.50" "1.60" "1.70" "1.80" "1.90" "2.00"];
%     filebase = "/data/scc/thobi/220201_ReducedSmapleStepDeltat/xy_s/sqrtN_";
    filebase = "/data/scc/thobi/201207_equilibration/xy_s/scale/sqrtN_";
end
linewidth=1.3;
T_vals=str2double(T_str);
N_T=numel(T_vals);

% filebase = "/data/scc/thobi/210715_LinearTimeSampling//mxy_3.00/sqrtN_";
filemid = "/T_";

for i_T = 1:N_T
    
    T = T_vals(i_T);
    fprintf("%.3f ",T)
    filenames(i_T) = filebase + sqrtN + filemid + T_str(i_T) + fileend;
    S{i_T}=load(filenames(i_T),"M_ang_MSD","averaging_times");
    filenames_shorttime(i_T) = filebase_shorttime + sqrtN + filemid + T_str(i_T) + fileend;
    S_shorttime{i_T}=load(filenames_shorttime(i_T),"M_ang_MSD","averaging_times");
end
fprintf("\n");

%% MSD Plot
c_map = linspecer(N_T);
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);

for i_T = 1:N_T
    
    T = T_vals(i_T);
    dispname = sprintf('$T = %.2f$',T);
    t=S{i_T}.averaging_times;
    MSD=S{i_T}.M_ang_MSD;
    t_st=S_shorttime{i_T}.averaging_times;
    MSD_st=S_shorttime{i_T}.M_ang_MSD;
    st_max=max(t_st);
    t_full=[t_st(t_st<=st_max),t(t>st_max)];
    MSD_full=[MSD_st(t_st<=st_max),MSD(t>st_max)];
    
    loglog(t_full,MSD_full,'DisplayName',dispname,...
        'Color',c_map(i_T,:),'LineWidth',linewidth);
    hold on;
    
end
xlim([1 Inf])
ylim([1e-8 Inf])
hlegend=legend('Location','SouthEast','Interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_annotation,...
    'NumColumns',2);
%     legpos=get(hlegend,"Position");
%     legpos=legpos+[0 .1 0 0];
%     set(hlegend,"Position",legpos)
hXLabel = xlabel('$t$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$\langle \phi_m^2(t)\rangle$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

%     text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');

figname=sprintf('%s/%s_angMSD_sqrtN_%d',pathbase,model,sqrtN);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');

