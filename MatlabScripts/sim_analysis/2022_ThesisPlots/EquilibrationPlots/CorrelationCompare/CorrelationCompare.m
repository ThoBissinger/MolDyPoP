%% initialization
clear
cd /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/EquilibrationPlots/CorrelationCompare
run('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/initialization_script.m')
pathbase=pwd;
model="mxy";
sqrtN = 256;

T_str = [".01" ".07" ".14" ".17" ".19"];
T_str = [".03" ".05" ".09" ".11" ".13"];
T_str = [".01" ".03" ".05" ".07" ".09" ".11" ".13" ".14" ".17" ".19" ".21" ".23"];
% T_str = [".01" ".03" ".05" ".07" ".13"];
% T_str = [".09" ".11" ".14" ".17" ".19" ".21"];
T_vals = str2double(T_str);
N_T = numel(T_vals);
filebase = "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_" + sqrtN + "/T_";
% filebase = "/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00/sqrtN_" + sqrtN + "/T_";
fileend = "/samp_Dynamics_collect.mat";
for i_T = 1:N_T
    filename = filebase + T_str(i_T) + fileend;
    S{i_T} = load(filename,"averaging_times","qbin","gmperpmperp","gmperpmperp_collect","SCF_Spin_av_collect");
end

%% Fitting (for FFT)
gamma_vals=zeros(N_T,6,4);
om_1_vals=zeros(N_T,6,4);
for i_T = 1:N_T
    T = T_vals(i_T);
    t = S{i_T}.averaging_times;
    
    q_vals = S{i_T}.qbin;
    n_t = numel(t);
    n_q = numel(q_vals);
    for i_integer_q = 1:6
        [~,i_q] = min(abs(q_vals/q_vals(1) - i_integer_q));
        q = q_vals(i_q);
        for k = 1:4
            ind_run=(k-1)*125+(1:125);
            cf = real(mean(S{i_T}.gmperpmperp_collect(ind_run,i_q:n_q:end)));
            cf = cf/cf(1);
            fitob=fit(t',cf',"exp(-gamma*x/2)*(cos(omega_1*x)+.5*gamma/omega_1 * sin(omega_1 * x))","StartPoint",[1e-2,1e-2],"Lower",[0,0]);
            gamma_vals(i_T,i_integer_q,k)=fitob.gamma;
            om_1_vals(i_T,i_integer_q,k)=fitob.omega_1;
        end
    end
end

%% Cmperp vs t for different run selections
c_map = linspecer(4);

for i_T = 1:N_T
    T = T_vals(i_T);
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 1.2*columnwidth_cm 1.44*columnwidth_cm]);
    t = S{i_T}.averaging_times;
    
    q_vals = S{i_T}.qbin;
    n_t = numel(t);
    n_q = numel(q_vals);
    for i_integer_q = 1:6
        [~,i_q] = min(abs(q_vals/q_vals(1) - i_integer_q));
        q = q_vals(i_q);
        for k = 1:4
            if i_integer_q == 1
                vis_handle = "on";
            else
                vis_handle = "off";
            end
            ind_run=(k-1)*125+(1:125);
            cf = real(mean(S{i_T}.gmperpmperp_collect(ind_run,i_q:n_q:end)));
            cf = cf/cf(1);
            
%             dispname = sprintf('run $%d - %d$',ind_run(1),ind_run(end));
            dispname = sprintf('$t_{\\textrm{i}} = %d \\cdot 10^4$',(k-1));
    %         dispname = sprintf('$\\tau_{\\textrm{anneal}} = %.2g$',step_vals(i_step));
            plot(t,cf + 3*i_integer_q,'DisplayName',dispname,...
                'Color',c_map(k,:),...
                'LineWidth',1,...
                "HandleVisibility",vis_handle);
            hold on;
%             htext{k} = text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');
        end
        htext{k} = text(8e3,3*i_integer_q - 1.5,sprintf("q = %.3f",q),'FontSize',fontsize_axis,'FontName', 'cmr12');
    end
    ylim([0 24]);
    hlegend=legend('Location','NorthEast','Interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_annotation,...
        'NumColumns',3);
%     legpos=get(hlegend,"Position");
%     legpos=legpos+[0 .1 0 0];
%     set(hlegend,"Position",legpos)
    hXLabel = xlabel('$t$','interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    hYLabel = ylabel('$C_{m\perp}(q,t)/\chi_{m\perp}(q)$','interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    set(gca,"ytick",[])

    pos=get(gca,'position');
    htitle=title(sprintf('$N = (%d)^2$, $T = %.2f$',sqrtN,T),'Interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_titles);
%     text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');
    
    figname=sprintf('%s/%s_sqrtN_%d_T_%s_cf_vs_t',pathbase,model,sqrtN,T_str(i_T));
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');

    htitle=title(sprintf('$N = (%d)^2$, $T = %.2f$, short time',sqrtN,T),'Interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_titles);
    xlim([0 max(t)/10]);
    figname=sprintf('%s/%s_sqrtN_%d_T_%s_cf_vs_t_zoomed',pathbase,model,sqrtN,T_str(i_T));
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end


%% Smperp vs t for different run selections
c_map = linspecer(4);

for i_T = 1:N_T
    T = T_vals(i_T);

    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 1.2*columnwidth_cm 1.44*columnwidth_cm]);
    t = S{i_T}.averaging_times;
    
    q_vals = S{i_T}.qbin;
    n_t = numel(t);
    n_q = numel(q_vals);
    for i_integer_q = 1:4
        [~,i_q] = min(abs(q_vals/q_vals(1) - i_integer_q));
        q = q_vals(i_q);
        ft_run_max = 0;

        damp_rate = min(4/max(gamma_vals(i_T,i_integer_q,:)),max(t)/2);
        res_vals = exp(-t.^2/2/damp_rate^2);
        for k = 1:4
            if i_integer_q == 1
                vis_handle = "on";
            else
                vis_handle = "off";
            end
            ind_run=(k-1)*125+(1:125);
            cf = real(mean(S{i_T}.gmperpmperp_collect(ind_run,i_q:n_q:end)));
            cf = cf/cf(1);
            
            [ft_vals,om_vals] = FT_correlation(t, cf .* res_vals, 1e6);
            ft_vals = real(ft_vals);
%             dispname = sprintf('run $%d - %d$',ind_run(1),ind_run(end));
            dispname = sprintf('$t_{\\textrm{i}} = %d \\cdot 10^4$',(k-1));
    %         dispname = sprintf('$\\tau_{\\textrm{anneal}} = %.2g$',step_vals(i_step));
            plot(om_vals,ft_vals,'DisplayName',dispname,...
                'Color',c_map(k,:),...
                'LineWidth',1,...
                "HandleVisibility",vis_handle);
            hold on;
            [maxval,imax_cur] = max(ft_vals);
            if maxval > ft_run_max
                ft_run_max=maxval;
                imax=imax_cur;
            end
%             htext{k} = text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');
        end
        
        maxes(i_integer_q) = ft_run_max;
        om_maxes(i_integer_q)=abs(om_vals(imax));
        if om_maxes(i_integer_q) < .1*q
            [~,i_hwhm]=min(abs(ft_vals - maxval/2));
            htext{k} = text(abs(om_vals(i_hwhm)),ft_run_max,...
                sprintf("q = %.3f",q),...
                'HorizontalAlignment','left',...
                'FontSize',fontsize_axis,'FontName', 'cmr12');
        else
            htext{k} = text(abs(om_vals(imax)),1.2*ft_run_max,...
                sprintf("q = %.3f",q),...
                'HorizontalAlignment','center',...
                'FontSize',fontsize_axis,'FontName', 'cmr12');
        end
    end
    om_max=max(om_maxes);
    if om_max < .1*q
        xlim([-.2 .2]);
    else
        xlim([0 1.3*om_max]);
    end
    
    ylim([0 1.5*max(maxes)]);
%     ylim([0 24]);
    hlegend=legend('Location','NorthEast','Interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_annotation,...
        'NumColumns',3);
%     legpos=get(hlegend,"Position");
%     legpos=legpos+[0 .1 0 0];
%     set(hlegend,"Position",legpos)
    hXLabel = xlabel('$\omega$','interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    hYLabel = ylabel('$S_{m\perp}(q,\omega)/\chi_{m\perp}(q)$','interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    htitle=title(sprintf('$N = (%d)^2$, $T = %.2f$',sqrtN,T),'Interpreter','latex', ...
        'FontName', 'cmr12','FontSize', fontsize_titles);
%     text(0.025,0.95,"a)",'Units','normalized','FontSize',fontsize_labels,'FontName', 'cmr12');
    
    figname=sprintf('%s/%s_sqrtN_%d_T_%s_FFT_vs_t',pathbase,model,sqrtN,T_str(i_T));
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');

end
