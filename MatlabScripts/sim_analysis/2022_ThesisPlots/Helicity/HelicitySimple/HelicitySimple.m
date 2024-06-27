%% initialization
clear
pathbase="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/Helicity/HelicitySimple";
cd(pathbase)
run('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/initialization_script.m')
% pathbase=pwd;
model="xy";


if model == "mxy"
    sqrtN_vals = [16,32,64,128,256];
    T_str = [".01" ".03" ".05" ".07" ".09" ".11" ".13" ".14" ".15" ".155" ".16" ".165" ".17" ".175" ".18" ".185" ".19" ".195" ".20" ".205" ".21" ".22" ".23" ".24" ".25" ".27" ".29" ".31" ".33" ".35" ".37" ".40" ".43" ".46" ".49" ".52" ];
    filebase = "/data/scc/thobi/210715_LinearTimeSampling//mxy_3.00/sqrtN_";
    xmax_stable = .19;
    xmax=.3;
    H_sign = -1;
    fileend = "/samp_Dynamics_helicity";
    Ups_prefac = 3;
    T_KT=.173;
    H_x_ylim=[.05 .14];
elseif model == "xy"
    sqrtN_vals = [16,32,64,128];
%     T_str = [".10" ".20" ".30" ".40" ".50" ".60" ".70" ".80" ".85" ".91" ".95" "1.00"];
    T_str = [".10" ".20" ".30" ".40" ".50" ".60" ".70" ".80" ".85" ".87" ".89" ".90" ".91" ".93" ".95" ".97" "1.00" "1.03" "1.06" "1.09" "1.10" "1.12" "1.15" "1.18" "1.20" "1.21" "1.24" "1.30" "1.40" "1.50" "1.60" "1.70" "1.80" "1.90" "2.00"];
%     filebase = "/data/scc/thobi/220201_ReducedSmapleStepDeltat/xy_s/sqrtN_";
    filebase = "/data/scc/thobi/201207_equilibration/xy_s/scale/sqrtN_";
    filebase_anneal = "/data/scc/thobi/201207_equilibration/xy_s/anneal/sqrtN_";
    xmax=1.5;
    xmax_stable = 1;
    H_sign = 1;
    fileend = "/samp_eq_helicity";
    Ups_prefac = 1;
    T_KT=.89;
    H_x_ylim=[.6 2];
end
linewidth=1.3;
T_vals=str2double(T_str);
N_N=numel(sqrtN_vals);
N_T=numel(T_vals);

% filebase = "/data/scc/thobi/210715_LinearTimeSampling//mxy_3.00/sqrtN_";
filemid = "/T_";

rho_s_vals=zeros(N_N,N_T);
H_x_vals=zeros(N_N,N_T);
H_y_vals=zeros(N_N,N_T);
I_x_vals=zeros(N_N,N_T);
I_y_vals=zeros(N_N,N_T);
H_r_vals=zeros(N_N,N_T);
H_s_vals=zeros(N_N,N_T);
I_x_2_vals=zeros(N_N,N_T);
I_y_2_vals=zeros(N_N,N_T);

for i_N = 1:N_N
    
    sqrtN=sqrtN_vals(i_N);
    fprintf("%d ",sqrtN)
    for i_T = 1:N_T
        T = T_vals(i_T);
        fprintf("%.3f ",T)
        if model == "xy" && sqrtN == 16
            filenames(i_N,i_T) = filebase_anneal + sqrtN + filemid + T_str(i_T) + fileend;
        else
            filenames(i_N,i_T) = filebase + sqrtN + filemid + T_str(i_T) + fileend;
        end
        S{i_N,i_T}=load(filenames(i_N,i_T));

        rho_s_vals(i_N,i_T)=Ups_prefac*S{i_N,i_T}.rho_s;
        H_x_vals(i_N,i_T)=H_sign*S{i_N,i_T}.H_x;
        H_y_vals(i_N,i_T)=H_sign*S{i_N,i_T}.H_y;
        I_x_vals(i_N,i_T)=S{i_N,i_T}.I_x;
        I_y_vals(i_N,i_T)=S{i_N,i_T}.I_y;
        H_r_vals(i_N,i_T)=S{i_N,i_T}.H_r;
        H_s_vals(i_N,i_T)=S{i_N,i_T}.H_s;
        I_x_2_vals(i_N,i_T)=S{i_N,i_T}.I_x_2;
        I_y_2_vals(i_N,i_T)=S{i_N,i_T}.I_y_2;

    end
    fprintf("\n");
end

%% helicity 
c_map = linspecer(N_N);
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);

for i_N = 1:N_N
    
    sqrtN = sqrtN_vals(i_N);
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    Ups_calc=Ups_prefac*1/sqrtN^2/2*((H_x_vals(i_N,:)+H_y_vals(i_N,:)) ...
        - 1./T_vals.*(I_x_2_vals(i_N,:)+I_y_2_vals(i_N,:)));
    plot(T_vals,Ups_calc,'DisplayName',dispname,...
        'Color',c_map(i_N,:),'LineWidth',linewidth);
    hold on;
    
end
xlim([0 xmax])
hlegend=legend('Location','SouthWest','Interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_annotation);
%     legpos=get(hlegend,"Position");
%     legpos=legpos+[0 .1 0 0];
%     set(hlegend,"Position",legpos)
hXLabel = xlabel('$T$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$\Upsilon$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

%     text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');

figname=sprintf('%s/%s_helicity',pathbase,model);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');


%% H_x 
c_map = linspecer(N_N);
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);

for i_N = 1:N_N
    
    sqrtN = sqrtN_vals(i_N);
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    plot(T_vals,(H_x_vals(i_N,:)+H_y_vals(i_N,:))/sqrtN^2,'DisplayName',dispname,...
        'Color',c_map(i_N,:),'LineWidth',linewidth);
    hold on;
    
end
xlim([0 xmax])
ylim(H_x_ylim);
hlegend=legend('Location','SouthWest','Interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_annotation);
%     legpos=get(hlegend,"Position");
%     legpos=legpos+[0 .1 0 0];
%     set(hlegend,"Position",legpos)
hXLabel = xlabel('$T$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$\langle H_x + H_y \rangle / 2N$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

%     text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');

figname=sprintf('%s/%s_H_x',pathbase,model);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');

%% I_x_2 
c_map = linspecer(N_N);
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);

for i_N = 1:N_N
    
    sqrtN = sqrtN_vals(i_N);
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    plot(T_vals,(I_x_2_vals(i_N,:)+I_y_2_vals(i_N,:))/sqrtN^2,'DisplayName',dispname,...
        'Color',c_map(i_N,:),'LineWidth',linewidth);
    hold on;
    
end
xlim([0 xmax])
hlegend=legend('Location','NorthWest','Interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_annotation);
%     legpos=get(hlegend,"Position");
%     legpos=legpos+[0 .1 0 0];
%     set(hlegend,"Position",legpos)
hXLabel = xlabel('$T$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$\langle I_x^2 + I_y^2 \rangle / 2N$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

%     text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');

figname=sprintf('%s/%s_I_x_2',pathbase,model);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');

%% beta I_x_2 
c_map = linspecer(N_N);
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);

for i_N = 1:N_N
    
    sqrtN = sqrtN_vals(i_N);
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    plot(T_vals,1./T_vals.*(I_x_2_vals(i_N,:)+I_y_2_vals(i_N,:))/sqrtN^2,'DisplayName',dispname,...
        'Color',c_map(i_N,:),'LineWidth',linewidth);
    hold on;
    
end
xlim([0 xmax])
hlegend=legend('Location','NorthWest','Interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_annotation);
%     legpos=get(hlegend,"Position");
%     legpos=legpos+[0 .1 0 0];
%     set(hlegend,"Position",legpos)
hXLabel = xlabel('$T$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$\beta \langle I_x^2 + I_y^2 \rangle / 2N$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

%     text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');

figname=sprintf('%s/%s_beta_I_x_2',pathbase,model);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');


%% helicity 
c_map = linspecer(N_N);
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);

for i_N = 1:N_N
    
    sqrtN = sqrtN_vals(i_N);
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    Ups_calc=Ups_prefac*1/sqrtN^2/2*((H_x_vals(i_N,:)+H_y_vals(i_N,:)) ...
        - 1./T_vals.*(I_x_2_vals(i_N,:)+I_y_2_vals(i_N,:)));
    plot(T_vals,T_vals./Ups_calc/2/pi,'DisplayName',dispname,...
        'Color',c_map(i_N,:),'LineWidth',linewidth);
    hold on;
    
end
xlim([0 xmax])
ylim([0 .4]);
xline(T_KT,'k--','HandleVisibility','off');
% ylim([0 .3])
hlegend=legend('Location','NorthWest','Interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_annotation);
%     legpos=get(hlegend,"Position");
%     legpos=legpos+[0 .1 0 0];
%     set(hlegend,"Position",legpos)
hXLabel = xlabel('$T$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$\eta$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);

text(.96*T_KT,0.03,"$T_{\textrm{BKT}}$",...
    'HorizontalAlignment','right','VerticalAlignment','bottom',...
    'interpreter','latex',...
    'FontSize',fontsize_annotation,'FontName', 'cmr12');

figname=sprintf('%s/%s_eta',pathbase,model);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');

xlim([0 xmax_stable])
figname=sprintf('%s/%s_eta_stable',pathbase,model);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');