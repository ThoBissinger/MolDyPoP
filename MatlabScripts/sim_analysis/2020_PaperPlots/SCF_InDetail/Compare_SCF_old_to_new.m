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
    T_str=[".16", ".17", ".18", ".19", ".20", ".21"];
%     T_str = [".03" ".05" ".07" ".09" ".11" ".13" ".14" ".15" ".155" ".16" ".165" ".17" ".175" ".18" ".185" ".19" ".195" ".20" ".205" ".21" ".22" ".23" ".24" ".25"];
    T_str=[".14"];
    T_vals=str2double(T_str);
    filedir="/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00";
    filename="samp_Dynamics_grSCF.mat";
    filename_old="samp_Dynamics.mat";
end
i_N=4;
sqrtN=sqrtN_vals(i_N);
L=L_vals(i_N);
n_T = numel(T_vals);
for i_T=1:n_T
    T=T_vals(i_T);
    curfile=sprintf("%s/sqrtN_%d/T_%s/%s",filedir,sqrtN,T_str(i_T),filename);
    curfile_old=sprintf("%s/sqrtN_%d/T_%s/%s",filedir,sqrtN,T_str(i_T),filename_old);
    S=load(curfile);
    S_old=load(curfile_old);

    %% Figure: C_m new and old
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    plot(S.rbin,S.Cm./S.gr,'r',...
        'DisplayName',"$C_{m}^{\textrm{new}}(r)$");
    hold on;
    plot(S_old.rbin,S_old.SCF_Spin_av/1.3,'bo',...
        'DisplayName',"$C_{m}^{\textrm{old}}(r)$");
    ylim([0 1]);
    xlim([1 Inf]);

    hXLabel = xlabel('$r$','interpreter','latex');
    hYLabel = ylabel('$C_{m}(r)$','interpreter','latex');
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
        'YScale','log', 'XScale','log',...
        'XAxisLocation','bottom',...
        'LineWidth', .5)
    


    figname=sprintf('%s/%s_cross_diff',basedir,model);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        %                 print(sprintf('%s.eps',figname),'-depsc');
    end



    %% Figure: C_m new by Cm_old
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    plot(S.rbin,S.Cm./S_old.SCF_Spin_av,'r',...
        'DisplayName',"$C_{m}^{\textrm{new}}(r)$");
    hold on;
    ylim([0 1]);
    xlim([1 Inf]);

    hXLabel = xlabel('$r$','interpreter','latex');
    hYLabel = ylabel('$C_{m}(r)$','interpreter','latex');
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
        'YScale','log', 'XScale','log',...
        'XAxisLocation','bottom',...
        'LineWidth', .5)
    


    figname=sprintf('%s/%s_cross_diff',basedir,model);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        %                 print(sprintf('%s.eps',figname),'-depsc');
    end
end