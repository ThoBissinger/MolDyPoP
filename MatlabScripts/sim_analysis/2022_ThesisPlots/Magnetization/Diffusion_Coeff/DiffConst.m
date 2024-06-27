%% initialization
clear
pathbase="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/Magnetization/Diffusion_Coeff";
cd(pathbase)
run('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/initialization_script.m')
% pathbase=pwd;
model="mxy";


if model == "mxy"
    sqrtN_vals = [16,32,64,128,256];
    T_str = [".01" ".03" ".05" ".07" ".09" ".11" ".13" ".14" ".15" ".155" ".16" ".165" ".17" ".175" ".18" ".185" ".19" ".195" ".20" ".205" ".21" ".22" ".23" ".24" ".25" ".27" ".29" ".31" ".33" ".35" ".37" ".40" ".43" ".46" ".49" ".52" ];
    filebase = "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_";
    fileend="/samp_Dynamics_M_TimeEvolution.mat";
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
N_N=numel(sqrtN_vals);
N_T=numel(T_vals);

% filebase = "/data/scc/thobi/210715_LinearTimeSampling//mxy_3.00/sqrtN_";
filemid = "/T_";

Diff_coeff=zeros(N_N,N_T);
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    fprintf("sqrtN = %d   ",sqrtN);
    for i_T = 1:N_T
        
        T = T_vals(i_T);
        fprintf("%.3f ",T)
        filenames(i_N,i_T) = filebase + sqrtN + filemid + T_str(i_T) + fileend;
        S{i_N,i_T}=load(filenames(i_N,i_T),"M_ang_MSD","averaging_times");
        Diff_coeff(i_N,i_T)=max(S{i_N,i_T}.M_ang_MSD)/max(S{i_N,i_T}.averaging_times);
    end
    fprintf("\n");
end

for i_T=1:N_T
    fitob{i_T}=fit(sqrtN_vals(:),Diff_coeff(:,i_T),"a*x^(-b)",...
        'StartPoint',[1,1],"Lower",[0 0]);
    Diff_a_vals(i_T)=fitob{i_T}.a;
    Diff_b_vals(i_T)=fitob{i_T}.b;
end
fprintf("\n");
savefile=sprintf("%s_data",model);


save(savefile)

%% MSD Plot
c_map = linspecer(N_N);
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);

for i_N = 1:N_N
    
    sqrtN = sqrtN_vals(i_N);
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    semilogy(T_vals,Diff_coeff(i_N,:),'DisplayName',dispname,...
        'Color',c_map(i_N,:),'LineWidth',linewidth);
    hold on;
    
end
xlim([0 .25]);
hlegend=legend('Location','SouthEast','Interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_annotation,...
    'NumColumns',2);
%     legpos=get(hlegend,"Position");
%     legpos=legpos+[0 .1 0 0];
%     set(hlegend,"Position",legpos)
hXLabel = xlabel('$T$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$D_m$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

% inset
innerpos=get(gca,'InnerPosition');
ax_inset = axes('Position',[innerpos(1)+.06 innerpos(2)+innerpos(4)*.6 innerpos(3)*.4 innerpos(4)*.38]);
plot(T_vals,Diff_b_vals,...
    'Color','r','LineWidth',linewidth);
xlim([0 .2]);
ylim([2 3])
text(.01,3,"$\sigma_D$","FontSize",fontsize_subfiglabels,"FontName","cmr12",...
    "Interpreter","latex",...
    "HorizontalAlignment","left","VerticalAlignment","top",...
    "Units","data");
text(.18,2.0,"$T$","FontSize",fontsize_subfiglabels,"FontName","cmr12",...
    "Interpreter","latex",...
    "HorizontalAlignment","right","VerticalAlignment","bottom",...
    "Units","data");

text(.2,3.0,"$D_m \sim L^{-\sigma_D}$","FontSize",fontsize_subfiglabels,"FontName","cmr12",...
    "Interpreter","latex",'Color','r',...
    "HorizontalAlignment","right","VerticalAlignment","top",...
    "Units","data");
%     text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');

figname=sprintf('%s/%s_DiffCoeff',pathbase,model);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');


%% Diff a/b
c_map = linspecer(N_N);
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);

plot(T_vals,Diff_b_vals,...
    'Color','k','LineWidth',linewidth);
    
xlim([0 .19]);
hXLabel = xlabel('$T$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$D_m$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

%     text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');

figname=sprintf('%s/%s_Diff_b',pathbase,model);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
