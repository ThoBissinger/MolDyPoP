close all;
clear all;


figbase='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/equilibration_checks';
saveswitch=1;

matfiles=["/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.05/samp_Dynamics.mat",...
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.09/samp_Dynamics.mat",...
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.11/samp_Dynamics.mat",...
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.155/samp_Dynamics.mat",...
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.17/samp_Dynamics.mat",...
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.19/samp_Dynamics.mat"];
T_str=[".05", ".09", ".11", ".155", ".17", ".19"];

runmax=500;
run_for_avg=125;

for j = 1:numel(matfiles)
    curmat=matfile(matfiles(j));
%     load(curmat,'averaging_times');
%     load(curmat,'ACF_Spin_collect');
    averaging_times=curmat.averaging_times;
    ACF_Spin_collect=curmat.ACF_Spin_collect;
    it_rounds=runmax/run_for_avg;
    
    c_map=linspecer(it_rounds);
    for i=1:it_rounds
        if (j == 1)
            visswitch = 'on';
        else
            visswitch = 'off';
        end
        av_ind=((i-1)*run_for_avg+1):(i*run_for_avg);
        dispname=sprintf('runs %d - %d',av_ind(1),av_ind(end));
%         plot(curmat.averaging_times,mean(curmat.ACF_Spin_collect(av_ind,:)),...
        h_plot{i,j}=plot(averaging_times,mean(ACF_Spin_collect(av_ind,:)),...
            'DisplayName',dispname,...
            'LineWidth',2,...
            'HandleVisibility',visswitch,...
            'Color',c_map(i,:)); 
        hold on; 
    end
    legend show
end
% legend([h_plot{1,1},h_plot{2,1},h_plot{3,1},h_plot{4,1}]);
ylim([0 1.2]);
xlabel('$t$','interpreter','latex');
ylabel('$C_m^{\textrm{ACF}}(t)$','interpreter','latex');

figname=sprintf('%s/mxy_equilibration_ACF_runcompare_longtime',figbase);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

xlim([0 2e3])
figname=sprintf('%s/mxy_equilibration_ACF_runcompare_shorttime',figbase);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end



figure(2)
c_map=linspecer(numel(matfiles));
for j = 1:numel(matfiles)
    curmat=matfile(matfiles(j));
%     load(curmat,'averaging_times');
%     load(curmat,'ACF_Spin_collect');
    averaging_times=curmat.averaging_times;
    ACF_Spin_collect=curmat.ACF_Spin_collect;
    it_rounds=runmax/run_for_avg;
    dispname=sprintf('T = %s',T_str(j));
    
    
    for i=1:it_rounds
        if (i == 1)
            visswitch = 'on';
        else
            visswitch = 'off';
        end
        av_ind=((i-1)*run_for_avg+1):(i*run_for_avg);
%         dispname=sprintf('runs %d - %d',av_ind(1),av_ind(end));
%         plot(curmat.averaging_times,mean(curmat.ACF_Spin_collect(av_ind,:)),...
        h_plot{i,j}=plot(averaging_times,mean(ACF_Spin_collect(av_ind,:)),...
            'DisplayName',dispname,...
            'LineWidth',2,...
            'HandleVisibility',visswitch,...
            'Color',c_map(j,:)); 
        hold on; 
    end
    legend show
end
% legend([h_plot{1,1},h_plot{1,2},h_plot{1,3},h_plot{1,4},h_plot{1,5},h_plot{1,6}]);
ylim([0 1.2]);
xlabel('$t$','interpreter','latex');
ylabel('$C_m^{\textrm{ACF}}(t)$','interpreter','latex');

figname=sprintf('%s/mxy_equilibration_ACF_Tcompare_longtime',figbase);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

xlim([0 2e3])
figname=sprintf('%s/mxy_equilibration_ACF_Tcompare_shorttime',figbase);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end