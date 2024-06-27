clear
curdir="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/SCF_InDetail";
cd(curdir);
run ../initialization_script
saveswitch=1;
basedir=sprintf('%s/Plots',curdir);

model = "mxy";
if model == "mxy"
%     sqrtN_vals=[16, 32, 64, 128, 256];
%     L_vals=9.25*2.^[0:4];
    sqrtN_vals=[16, 32, 64, 128];
    L_vals=9.25*2.^[0:3];
    T_str=[".09", ".17", ".20"];
%     T_str = [".03" ".05" ".07" ".09" ".11" ".13" ".14" ".15" ".155" ".16" ".165" ".17" ".175" ".18" ".185" ".19" ".195" ".20" ".205" ".21" ".22" ".23" ".24" ".25"];
    T_str=[".14"];
    T_vals=str2double(T_str);
    filedir="/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00";
    filename="samp_Dynamics_grSCF.mat";
end
i_N=4;
sqrtN=sqrtN_vals(i_N);
L=L_vals(i_N);
n_T = numel(T_vals);
for i_T=1:n_T
    T=T_vals(i_T);
    curfile=sprintf("%s/sqrtN_%d/T_%s/%s",filedir,sqrtN,T_str(i_T),filename);
    S=load(curfile);

    %% Figure: <m^2gr> and <s.mgr> by <m^2><gr>
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    plot(S.rbin,smoothen(S.Cm_cross/S.M_2-S.gr,0,0),'r',...
        'DisplayName',"$C_{m}^{\textrm{cross}}(r)$");
    hold on;
    plot(S.rbin,smoothen(S.gr_absM_2/S.M_2-S.gr,0,0),'b',...
        'DisplayName',"$C_{m}^{m^2}(r)$");
    xlim([0 1.5]);

    hXLabel = xlabel('$r$','interpreter','latex');
    hYLabel = ylabel('$\rho g(r)\left[\frac{C_{m}^{\alpha}(r)}{\langle m^2\rangle} - 1\right]$','interpreter','latex');
    hLegend = legend('Location', 'SouthEast','interpreter','latex',...
        'NumColumns',1);
    htitle=title(sprintf('$N=(%d)^2$, $T=%.3f$',sqrtN,T),'Interpreter','latex',...
        'FontName','cmr14','fontsize',fontsize_titles);
    
    h_axis = gca;
    % Font
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
    
    % Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'XGrid', 'off', 'YGrid', 'off', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'YScale','lin', 'XScale','lin',...
        'XAxisLocation','bottom',...
        'LineWidth', .5)
    


    figname=sprintf('%s/%s_cross_diff',basedir,model);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        %                 print(sprintf('%s.eps',figname),'-depsc');
    end

    %% Figure: <m^2gr> and <s.mgr> alongside <m^2><gr>, zoomed
    xmin=.3;
    xmax=1.2;
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    plot(S.rbin,smoothen(S.Cm_cross/S.M_2,0,0),'r',...
        'DisplayName',"$\rho g(r)C_{m}^{\textrm{cross}}(r)/\langle m^2\rangle$",...
        'LineWidth',.8);
    hold on;
    plot(S.rbin,smoothen(S.gr_absM_2/S.M_2,0,0),'b',...
        'DisplayName',"$\rho g(r)C_{m}^{m^2}(r)/\langle m^2\rangle$",...
        'LineWidth',.8);
    plot(S.rbin,smoothen(S.gr,0,0),'k',...
        'DisplayName',"$\rho g(r)$",...
        'LineWidth',.8);


    hXLabel = xlabel('$r$','interpreter','latex');
    hYLabel = ylabel('Correlation','interpreter','latex');
    hLegend = legend('Location', 'SouthEast','interpreter','latex',...
        'NumColumns',1);
    htitle=title(sprintf('$N=(%d)^2$, $T=%.3f$',sqrtN,T),'Interpreter','latex',...
        'FontName','cmr14','fontsize',fontsize_titles);
    
    h_axis = gca;
    % Font
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
    
    % Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'XGrid', 'off', 'YGrid', 'off', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'YScale','lin', 'XScale','lin',...
        'XAxisLocation','bottom',...
        'LineWidth', .5)
    xlim([xmin xmax]);
    ylim([.7 1.2]);

    figname=sprintf('%s/%s_cross_zoom',basedir,model);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        %                 print(sprintf('%s.eps',figname),'-depsc');
    end


    %% Figure: C_m vs C_m^0
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    plot(S.rbin(S.gr>1e-2),S.Cm(S.gr>1e-2)./S.gr(S.gr>1e-2),'r',...
        'DisplayName',"$C_{m}(r)$",...
        'LineWidth',.8);
    hold on;
    plot(S.rbin(S.gr>1e-2),(S.Cm(S.gr>1e-2)-2*S.Cm_cross(S.gr>1e-2)+S.gr_absM_2(S.gr>1e-2))./S.gr(S.gr>1e-2),'b',...
        'DisplayName',"$C_{m}^{0}(r)$",...
        'LineWidth',.8);


    hXLabel = xlabel('$r$','interpreter','latex');
    hYLabel = ylabel('Correlation','interpreter','latex');
    hLegend = legend('Location', 'East','interpreter','latex',...
        'NumColumns',1);
    htitle=title(sprintf('$N=(%d)^2$, $T=%.3f$',sqrtN,T),'Interpreter','latex',...
        'FontName','cmr14','fontsize',fontsize_titles);
    
    h_axis = gca;
    % Font
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
    
    % Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'XGrid', 'off', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'YScale','lin', 'XScale','lin',...
        'XAxisLocation','bottom',...
        'LineWidth', .5)
    
    figname=sprintf('%s/%s_Cm_vs_Cm0',basedir,model);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        %                 print(sprintf('%s.eps',figname),'-depsc');
    end

    %% Figure: C_m vs C_m_par
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    plot(S.rbin(S.gr>1e-2),S.Cm(S.gr>1e-2)./S.gr(S.gr>1e-2),'r',...
        'DisplayName',"$C_{m}(r)$",...
        'LineWidth',.8);
    hold on;
    plot(S.rbin(S.gr>1e-2),(S.Cm_par(S.gr>1e-2))./S.gr(S.gr>1e-2),'b',...
        'DisplayName',"$C_{m\parallel}(r)$",...
        'LineWidth',.8);


    hXLabel = xlabel('$r$','interpreter','latex');
    hYLabel = ylabel('Correlation','interpreter','latex');
    hLegend = legend('Location', 'East','interpreter','latex',...
        'NumColumns',1);
    htitle=title(sprintf('$N=(%d)^2$, $T=%.3f$',sqrtN,T),'Interpreter','latex',...
        'FontName','cmr14','fontsize',fontsize_titles);
    
    h_axis = gca;
    % Font
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
    
    % Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'XGrid', 'off', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'YScale','lin', 'XScale','lin',...
        'XAxisLocation','bottom',...
        'LineWidth', .5)
    
    figname=sprintf('%s/%s_Cm_vs_Cmpar',basedir,model);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        %                 print(sprintf('%s.eps',figname),'-depsc');
    end


    %% Figure: C_m_par^0 vs C_m_perp
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    plot(S.rbin(S.gr>1e-2),S.Cm_perp(S.gr>1e-2)./S.gr(S.gr>1e-2),'r',...
        'DisplayName',"$C_{m\perp}(r)$",...
        'LineWidth',.8);
    hold on;
    plot(S.rbin(S.gr>1e-2),(S.Cm_par(S.gr>1e-2)-2*S.Cm_cross(S.gr>1e-2)+S.gr_absM_2(S.gr>1e-2))./S.gr(S.gr>1e-2),'b',...
        'DisplayName',"$C_{m\parallel}^0(r)$",...
        'LineWidth',.8);


    hXLabel = xlabel('$r$','interpreter','latex');
    hYLabel = ylabel('Correlation','interpreter','latex');
    hLegend = legend('Location', 'East','interpreter','latex',...
        'NumColumns',1);
    htitle=title(sprintf('$N=(%d)^2$, $T=%.3f$',sqrtN,T),'Interpreter','latex',...
        'FontName','cmr14','fontsize',fontsize_titles);
    
    h_axis = gca;
    % Font
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
    
    % Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'XGrid', 'off', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'YScale','lin', 'XScale','lin',...
        'XAxisLocation','bottom',...
        'LineWidth', .5)
    
    figname=sprintf('%s/%s_Cm_vs_Cmpar',basedir,model);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        %                 print(sprintf('%s.eps',figname),'-depsc');
    end
end

%% A/eta/xi values (fits)
T_str = [".01" ".03" ".05" ".07" ".09" ".11" ".13" ".14" ".15" ".155" ".16" ".165" ".17" ".175" ".18" ".185" ".19" ".195" ".20" ".205" ".21" ".22" ".23" ".24" ".25" ".27" ".29" ".31" ".33" ".35" ".37" ".40" ".43" ".46" ".49" ".52"];
fitob_cell_pow_Cm=cell(numel(sqrtN_vals),numel(T_str));
fitob_cell_exp_Cm=cell(numel(sqrtN_vals),numel(T_str));
fitob_cell_exp_Cm_par=cell(numel(sqrtN_vals),numel(T_str));

Cm_end_vals=zeros(numel(sqrtN_vals),numel(T_str));
Cm_par_end_vals=zeros(numel(sqrtN_vals),numel(T_str));
eta_vals=zeros(numel(sqrtN_vals),numel(T_str));
a_vals=zeros(numel(sqrtN_vals),numel(T_str));
a_vals_at_one=zeros(numel(sqrtN_vals),numel(T_str));
absM_vals=zeros(numel(sqrtN_vals),numel(T_str));
M_2_vals=zeros(numel(sqrtN_vals),numel(T_str));
M_4_vals=zeros(numel(sqrtN_vals),numel(T_str));

xi_exp_Cm_vals=zeros(numel(sqrtN_vals),numel(T_str));
a_exp_Cm_vals=zeros(numel(sqrtN_vals),numel(T_str));
c_exp_Cm_vals=zeros(numel(sqrtN_vals),numel(T_str));

xi_exp_Cm_par_vals=zeros(numel(sqrtN_vals),numel(T_str));
a_exp_Cm_par_vals=zeros(numel(sqrtN_vals),numel(T_str));
c_exp_Cm_par_vals=zeros(numel(sqrtN_vals),numel(T_str));

for i_N=1:numel(sqrtN_vals)
    for i_T=1:numel(T_str)
        curfile=sprintf("%s/sqrtN_%d/T_%s/%s",filedir,sqrtN_vals(i_N),T_str(i_T),filename);
        S=load(curfile,'rbin','Cm','Cm_par','gr','absM','M_2','M_4');
        fitob_cell_pow_Cm{i_N,i_T}=fit(S.rbin(S.rbin>3)',S.Cm(S.rbin>3)'./S.gr(S.rbin>3)',"a*x^(-b)","StartPoint",[1,.1],"Lower",[0,0]);
        [~,imin]=min(abs(S.rbin-1));
        a_vals(i_N,i_T)=fitob_cell_pow_Cm{i_N,i_T}.a;
        a_vals_at_one(i_N,i_T)=S.Cm(imin)/S.gr(imin);
        eta_vals(i_N,i_T)=fitob_cell_pow_Cm{i_N,i_T}.b;
        Cm_end_vals(i_N,i_T)=S.Cm(end-1)/S.gr(end-1);
        Cm_par_end_vals(i_N,i_T)=S.Cm_par(end-1)/S.gr(end-1);
        absM_vals(i_N,i_T)=S.absM;
        M_2_vals(i_N,i_T)=S.M_2;
        M_4_vals(i_N,i_T)=S.M_4;

        fitob_cell_exp_Cm{i_N,i_T}=fit(S.rbin(S.rbin>3)',S.Cm(S.rbin>3)'./S.gr(S.rbin>3)',"a*exp(-x/b)+c","StartPoint",[1,1,0],"Lower",[0,0,-1],"Upper",[1,L_vals(i_N),1]);
        a_exp_Cm_vals(i_N,i_T)=fitob_cell_exp_Cm{i_N,i_T}.a;
        xi_exp_Cm_vals(i_N,i_T)=fitob_cell_exp_Cm{i_N,i_T}.b;
        c_exp_Cm_vals(i_N,i_T)=fitob_cell_exp_Cm{i_N,i_T}.c;

        fitob_cell_exp_Cm_par{i_N,i_T}=fit(S.rbin(S.rbin>3)',S.Cm_par(S.rbin>3)'./S.gr(S.rbin>3)',"a*exp(-x/b)+c","StartPoint",[1,1,0],"Lower",[0,0,-1],"Upper",[1,L_vals(i_N),1]);
        a_exp_Cm_par_vals(i_N,i_T)=fitob_cell_exp_Cm_par{i_N,i_T}.a;
        xi_exp_Cm_par_vals(i_N,i_T)=fitob_cell_exp_Cm_par{i_N,i_T}.b;
        c_exp_Cm_par_vals(i_N,i_T)=fitob_cell_exp_Cm_par{i_N,i_T}.c;
    end
end
%% eta values plot
close all
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
c_map=lines(numel(sqrtN_vals));
for i_N=1:numel(sqrtN_vals)
    dispname=sprintf("$N = (%d)^2$",sqrtN_vals(i_N));
    plot(str2double(T_str),eta_vals(i_N,:),'s-',...
        "DisplayName",dispname,...
        "Color",c_map(i_N,:))
    hold on;
end
xlim([0 .25]); ylim([0 .35]);
hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$\eta$','interpreter','latex');
hLegend = legend('Location', 'NorthWest','interpreter','latex',...
    'NumColumns',1);

h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','lin',...
    'XAxisLocation','bottom',...
    'LineWidth', .5)

figname=sprintf('%s/%s_eta',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end




%% a values plot
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
c_map=lines(numel(sqrtN_vals));
for i_N=1:numel(sqrtN_vals)
    dispname=sprintf("$N = (%d)^2$",sqrtN_vals(i_N));
    plot(str2double(T_str),a_vals(i_N,:),'s-',...
        "DisplayName",dispname,...
        "Color",c_map(i_N,:))
    hold on;
end
ylim([0 1]);
xlim([0 .5]); 

hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$A$','interpreter','latex');
hLegend = legend('Location', 'SouthWest','interpreter','latex',...
    'NumColumns',1);

h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','lin',...
    'XAxisLocation','bottom',...
    'LineWidth', .5)

figname=sprintf('%s/%s_a',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end


%% a by m values plot
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
c_map=lines(numel(sqrtN_vals));
for i_N=1:numel(sqrtN_vals)
    dispname=sprintf("$N = (%d)^2$",sqrtN_vals(i_N));
    plot(str2double(T_str),a_vals(i_N,:)./absM_vals(i_N,:),'s-',...
        "DisplayName",dispname,...
        "Color",c_map(i_N,:))
    hold on;
end
ylim([0 1.2]);
xlim([0 .25]); 

hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$A/\langle m\rangle$','interpreter','latex');
hLegend = legend('Location', 'SouthWest','interpreter','latex',...
    'NumColumns',1);

h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','lin',...
    'XAxisLocation','bottom',...
    'LineWidth', .5)

figname=sprintf('%s/%s_a_by_m',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end


%% Cm(1) by a values plot
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
c_map=lines(numel(sqrtN_vals));
for i_N=1:numel(sqrtN_vals)
    dispname=sprintf("$N = (%d)^2$",sqrtN_vals(i_N));
    plot(str2double(T_str),a_vals_at_one(i_N,:)./a_vals(i_N,:),'s-',...
        "DisplayName",dispname,...
        "Color",c_map(i_N,:))
    hold on;
end
ylim([0 1.2]);
xlim([0 .25]); 

hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$C_m(1)/A$','interpreter','latex');
hLegend = legend('Location', 'SouthWest','interpreter','latex',...
    'NumColumns',1);


h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','lin',...
    'XAxisLocation','bottom',...
    'LineWidth', .5)

figname=sprintf('%s/%s_Cm1_by_a',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end


%% a^(-eta) values plot
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
c_map=lines(numel(sqrtN_vals));
for i_N=1:numel(sqrtN_vals)
    dispname=sprintf("$N = (%d)^2$",sqrtN_vals(i_N));
    plot(str2double(T_str),a_vals(i_N,:).^(1./eta_vals(i_N,:)),'s-',...
        "DisplayName",dispname,...
        "Color",c_map(i_N,:))
    hold on;
end
ylim([0 4e-4]);
xlim([0 .25]); 

hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$A^{-\eta}$','interpreter','latex');
hLegend = legend('Location', 'NorthEast','interpreter','latex',...
    'NumColumns',1);


h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','lin',...
    'XAxisLocation','bottom',...
    'LineWidth', .5)

figname=sprintf('%s/%s_a_to_1_by_eta',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end


%% xi_parallel values plot
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
c_map=lines(numel(sqrtN_vals));
for i_N=1:numel(sqrtN_vals)
    dispname=sprintf("$N = (%d)^2$",sqrtN_vals(i_N));
    plot(str2double(T_str),xi_exp_Cm_par_vals(i_N,:),'s-',...
        "DisplayName",dispname,...
        "Color",c_map(i_N,:))
    hold on;
end
xlim([0 .5]); 
ylim([0 20]);

hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$\xi_{\parallel}$','interpreter','latex');
hLegend = legend('Location', 'NorthEast','interpreter','latex',...
    'NumColumns',1);


h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','lin',...
    'XAxisLocation','bottom',...
    'LineWidth', .5)

figname=sprintf('%s/%s_xi_par_full',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end

xlim([.12 .35])
ylim([0 10])
hchild=get(gca,"Children");
delete(hchild(end));
figname=sprintf('%s/%s_xi_par_zoomed',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end

%% xi values plot
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
c_map=lines(numel(sqrtN_vals));
for i_N=1:numel(sqrtN_vals)
    dispname=sprintf("$N = (%d)^2$",sqrtN_vals(i_N));
    plot(str2double(T_str),xi_exp_Cm_vals(i_N,:),'s-',...
        "DisplayName",dispname,...
        "Color",c_map(i_N,:))
    hold on;
end
xlim([0 .5]); 
ylim([0 20]);

hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$\xi$','interpreter','latex');
hLegend = legend('Location', 'NorthEast','interpreter','latex',...
    'NumColumns',1);


h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','lin',...
    'XAxisLocation','bottom',...
    'LineWidth', .5)

figname=sprintf('%s/%s_xi',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end


xlim([.15 .35])
ylim([0 10])
figname=sprintf('%s/%s_xi_zoomed',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end

%%


figure
plot(str2double(T_str),xi_exp_Cm_par_vals,'s-');
%     'DisplayName',sprintfc("$sqrtN=(%d)^2$",sqrtN_vals));
hold on;
plot(str2double(T_str),xi_exp_Cm_vals,'o-');
ylim([0 max(L)])

figure
plot(str2double(T_str),xi_exp_Cm_vals,'s-');

