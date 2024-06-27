%% initialization
clear
cd /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/EquilibrationPlots/RateDependenceAnneal
run('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/initialization_script.m')
pathbase=pwd;
model="mxy";
sqrtN_vals = [16,32,64,128];
% sqrtN_vals = [64];
step_str=[".01" ".02" ".05" ".1" ".2" ".5" "1" "2" "4" "8" "16"];
step_str=[".05" ".1" ".2" ".5" "1" "2" "4" "8" "16"];
% step_str=[".05" ".1" ".2" ".5" "1" "2" "4" "8"];
% step_str=["2" "4" "8"];
step_vals=str2double(step_str);
N_N=numel(sqrtN_vals);
N_step=numel(step_vals);

filebase = "/data/scc/thobi/220829_EquilibrationCheck/anneal/ann_step_";
filemid = "/aligned/mxy_3.00/sqrtN_";
fileend = "/T_.01/samp_eq";
for i_N = 1:N_N
    sqrtN=sqrtN_vals(i_N);
    for i_step = 1:N_step
        step = step_vals(i_step);
        filenames(i_N,i_step) = filebase + step_str(i_step) + filemid + sqrtN + fileend;
        if step == 2
            filenames(i_N,i_step) = "/data/scc/thobi/220829_EquilibrationCheck/anneal" + filemid + sqrtN + fileend;
        end
        S{i_N,i_step}=load(filenames(i_N,i_step));
    end
end

%% m vs t dependence
c_map = turbo(N_step);

for i_N = 1:N_N
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);

    sqrtN = sqrtN_vals(i_N);
    for i_step = 1:N_step
        
        t = S{i_N,i_step}.averaging_times;
        m = S{i_N,i_step}.absM;
        dispname = sprintf('$\\tau = %.2g$',step_vals(i_step));
%         dispname = sprintf('$\\tau_{\\textrm{anneal}} = %.2g$',step_vals(i_step));
        semilogx(t,m,'DisplayName',dispname,...
            'Color',c_map(i_step,:));
        hold on;
    end
    xlim([0 1e7]);
    ylim([ 0 1]);
    hlegend=legend('Location','SouthEast','Interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_annotation);
%     legpos=get(hlegend,"Position");
%     legpos=legpos+[0 .1 0 0];
%     set(hlegend,"Position",legpos)
    hXLabel = xlabel('$t$','interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    hYLabel = ylabel('$\langle m \rangle$','interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    
%     text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');
    
    figname=sprintf('%s/%s_sqrtN_%d_m_vs_time',pathbase,model,sqrtN);
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

%% m vs temperature dependence
c_map = turbo(N_step);

for i_N = 1:N_N
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.2*columnwidth_cm]);

    sqrtN = sqrtN_vals(i_N);
    for i_step = 1:N_step
        
        T = S{i_N,i_step}.temperature(2:end);
        m = S{i_N,i_step}.absM(2:end);
%         T = smoothen(T,2,1000);
%         m = smoothen(m,2,1000);
        dispname = sprintf('$\\tau = %.2g$',step_vals(i_step));
%         dispname = sprintf('$\\tau_{\\textrm{anneal}} = %.2g$',step_vals(i_step));
        plot(T,m,'DisplayName',dispname,...
            'Color',c_map(i_step,:),...
            'LineWidth',1.5);
        hold on;
    end
    xlim([0 .5]);
    ylim([ 0 1]);
    hlegend=legend('Location','NorthEast','Interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_annotation);
    htitle=title(sprintf('$N = (%d)^2$',sqrtN),'Interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_titles);
%     legpos=get(hlegend,"Position");
%     legpos=legpos+[0 .1 0 0];
%     set(hlegend,"Position",legpos)
    hXLabel = xlabel('$\langle T \rangle$','interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    hYLabel = ylabel('$\langle m \rangle$','interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    
%     text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');
    
    figname=sprintf('%s/%s_sqrtN_%d_m_vs_temp',pathbase,model,sqrtN);
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

