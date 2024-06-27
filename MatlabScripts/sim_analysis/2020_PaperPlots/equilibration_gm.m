close all;
clear all;


figbase='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/equilibration_checks';
saveswitch=1;
mode="runcompare"; % Different runs have different colors
% mode="qcompare"; % Different wavevectors have different colors

% matfiles=["/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.01/samp_Dynamics.mat",...
%     "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.03/samp_Dynamics.mat",...
%     "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.05/samp_Dynamics.mat",...
%     "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.07/samp_Dynamics.mat",...
%     "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.09/samp_Dynamics.mat",...
%     "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.11/samp_Dynamics.mat",...
%     "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.155/samp_Dynamics.mat",...
%     "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.16/samp_Dynamics.mat",...
%     "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.165/samp_Dynamics.mat",...
%     "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.17/samp_Dynamics.mat",...
%     "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.175/samp_Dynamics.mat",...
%     "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.19/samp_Dynamics.mat"];
T_str=[".01", ".03", ".05", ".07", ".09", ".11", ".155", ".16", ".165", ".17", ".175", ".18", ".185", ".19"];
N_T=numel(T_str);
for i_T=1:N_T
    matfiles(i_T)=sprintf("/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_%s/samp_Dynamics.mat",T_str(i_T));
end
runmax=500;
run_for_avg=125;
q_select=1:6;
N_q=numel(q_select);
it_rounds=runmax/run_for_avg;
if ( mode == "qcompare") 
    c_map=linspecer(N_q);
elseif ( mode == "runcompare" )
    c_map=linspecer(it_rounds);
end
% for i_T = 1:N_T
for i_T = [12,13]
    figure(i_T)
    curmat=matfile(matfiles(i_T));
%     load(curmat,'averaging_times');
%     load(curmat,'ACF_Spin_collect');
    averaging_times=curmat.averaging_times;
    q_vals=curmat.qbin;
    gmperpmperp_collect=curmat.gmperpmperp_collect;

    
    for i_q = 1:N_q
        q_indices = (i_q):length(q_vals):numel(gmperpmperp_collect(1,:));
        q=q_vals(i_q);
        for i_run=1:it_rounds
            av_ind=((i_run-1)*run_for_avg+1):(i_run*run_for_avg);
            if ( mode == "qcompare" && i_run == 1) || ( mode == "runcompare" && i_q == 1)
                visswitch = 'on';
                
            else
                visswitch = 'off';
            end
            if ( mode == "qcompare") 
                dispname=sprintf('q = %.3f',q);
                color=c_map(i_q,:);
            elseif ( mode == "runcompare" )
                dispname=sprintf('runs %d - %d',av_ind(1),av_ind(end));
                color=c_map(i_run,:);
            end

%             dispname=sprintf('q = %.3f',q);
%             dispname=sprintf('runs %d - %d',av_ind(1),av_ind(end));
            h_plot{i_run,i_q}=plot(averaging_times,mean(real(gmperpmperp_collect(av_ind,q_indices)))/mean(real(gmperpmperp_collect(av_ind,i_q))) + 3*i_q,...
                'DisplayName',dispname,...
                'LineWidth',2,...
                'HandleVisibility',visswitch,...
                'Color',color);
            hold on;
        end
    end
    legend show
    set(legend,'NumColumns',3)
% legend([h_plot{1,1},h_plot{2,1},h_plot{3,1},h_plot{4,1}]);
    ylim([0 24]);
    xlabel('$t$','interpreter','latex');
    ylabel('$C_{m\perp}(q,t)/\chi_{m\perp}$','interpreter','latex');
    
    pos=get(gca,'position');
    annotation('textbox',[pos(1)+.04*pos(3),pos(2)+.845*pos(4),.1,.1],'FitBoxToText','on',...
        'string',sprintf('$T = %s$', T_str(i_T)),'interpreter','latex');

    figname=sprintf('%s/mxy_equilibration_gmperp_T_%s_%s_longtime',figbase,T_str(i_T),mode);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    xlim([0 2e3])
    figname=sprintf('%s/mxy_equilibration_gmperp_T_%s_%s_shorttime',figbase,T_str(i_T),mode);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end

return
figure(2)
c_map=linspecer(numel(matfiles));
for i_T = 1:numel(matfiles)
    curmat=matfile(matfiles(i_T));
%     load(curmat,'averaging_times');
%     load(curmat,'ACF_Spin_collect');
    averaging_times=curmat.averaging_times;
    ACF_Spin_collect=curmat.ACF_Spin_collect;
    it_rounds=runmax/run_for_avg;
    dispname=sprintf('T = %s',T_str(i_T));


    for i_run=1:it_rounds
        if (i_run == 1)
            visswitch = 'on';
        else
            visswitch = 'off';
        end
        av_ind=((i_run-1)*run_for_avg+1):(i_run*run_for_avg);
%         dispname=sprintf('runs %d - %d',av_ind(1),av_ind(end));
%         plot(curmat.averaging_times,mean(curmat.ACF_Spin_collect(av_ind,:)),...
        h_plot{i_run,i_T}=plot(averaging_times,mean(ACF_Spin_collect(av_ind,:)),...
            'DisplayName',dispname,...
            'LineWidth',2,...
            'HandleVisibility',visswitch,...
            'Color',c_map(i_T,:));
        hold on;
    end
    legend show

% legend([h_plot{1,1},h_plot{1,2},h_plot{1,3},h_plot{1,4},h_plot{1,5},h_plot{1,6}]);
% ylim([0 1.2]);
xlabel('$t$','interpreter','latex');
ylabel('$C_m(q,t)$','interpreter','latex');

figname=sprintf('%s/mxy_equilibration_Tcompare_longtime',figbase);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

xlim([0 2e3])
figname=sprintf('%s/mxy_equilibration_Tcompare_shorttime',figbase);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end
end