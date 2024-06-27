%% initialization
clear
pathbase="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/Magnetization/Susceptibility";
cd(pathbase)
run('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/initialization_script.m')
% pathbase=pwd;
model="mxy";


if model == "mxy"
    sqrtN_vals = [16,32,64,128,256];
    T_str = [".01" ".03" ".05" ".07" ".09" ".11" ".13" ".14" ".15" ".155" ".16" ".165" ".17" ".175" ".18" ".185" ".19" ".195" ".20" ".205" ".21" ".22" ".23" ".24" ".25" ".27" ".29" ".31" ".33" ".35" ".37" ".40" ".43" ".46" ".49" ".52" ];
    filebase = "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_";
    fileend="/samp_Dynamics.mat";
    filename="samp_Dynamics";
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

m_vals=zeros(N_N,N_T);
m_2_vals=zeros(N_N,N_T);
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    fprintf("sqrtN = %d   ",sqrtN);
    for i_T = 1:N_T
        
        T = T_vals(i_T);
        fprintf("%.3f ",T)
        filenames(i_N,i_T) = filebase + sqrtN + filemid + T_str(i_T) + fileend;
        S{i_N,i_T}=load(filenames(i_N,i_T),"absM_av","M_2_av");
        m_vals(i_N,i_T)=S{i_N,i_T}.absM_av;
        m_2_vals(i_N,i_T)=S{i_N,i_T}.M_2_av;
    end
    fprintf("\n");
end
sigma_vals=(m_2_vals - m_vals.^2);
chi_vals=sqrtN_vals'.^2./T_vals.*sigma_vals;


for i_T=1:N_T
    fitob{i_T}=fit(sqrtN_vals(:),chi_vals(:,i_T)./sqrtN_vals',"a*x^(-b)",...
        'StartPoint',[1,1],"Lower",[0 0]);
    a_vals(i_T)=fitob{i_T}.a;
    eta_vals(i_T)=fitob{i_T}.b;
end
fprintf("\n");
savefile=sprintf("%s_data",model);


save(savefile)



%% Susceptibility Plot
c_map = linspecer(N_N);
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);

for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    N = sqrtN^2;
    dispname=sprintf("N = (%d)^2",sqrtN);
    plot(T_vals,chi_vals(i_N,:),'--s',...
        "DisplayName",dispname,...
        "LineWidth",1.2,...
        'MarkerSize',4,...
        'Color',c_map(i_N,:),'LineWidth',linewidth);
    hold on;
end
xlim([0 .36]);

% hlegend = legend('Location',',...
%     'interpreter','latex', ...
%     'FontName', 'cmr12','FontSize', fontsize_legend);

hXLabel = xlabel('$T$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$\chi_m$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

%     text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');

figname=sprintf('%s/%s_Susceptibility',pathbase,model);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
