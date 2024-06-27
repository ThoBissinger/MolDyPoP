%% Initialization
clear all
% addpath "/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/";
run ../initialization_script;
pathbase="/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00";
figbase="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/BondOrder";
cd(figbase);
sqrtN=32;
L=9.25*sqrtN/16;
rho=(sqrtN/L)^2;
T_str=[".01" ".03" ".05" ".07" ".09" ".11" ".13" ".15"];
% T_str=[ ".01"  ".03"  ".09"  ".15"  ".19"  ".31"];
T_vals=str2double(T_str);
N_T=numel(T_vals);
rmax_vals=[0 0.66 0.73 0.80 1.08];
N_rmax=numel(rmax_vals);
pathbase=pathbase + "/sqrtN_" + sqrtN;
i_r=2;
r_trig=sqrt(2/sqrt(3)/rho);
r_hon=sqrt(2)*r_trig;

S=cell(N_T,N_rmax);
for i_T = 1:N_T
    curfile=pathbase+"/T_" + T_str(i_T) + "/samp_Dynamics_bondorder.mat";
    S_cur=load(curfile,"Psibin","Psihist_mean","rbin","Gcorrs_mean","rmax_vals","gr");
    for i_r = 1:N_rmax
        S{i_T,i_r}.Psibin=S_cur.Psibin;
        S{i_T,i_r}.Psi=reshape(S_cur.Psihist_mean(4,i_r,:),1,[]);
        S{i_T,i_r}.rbin=S_cur.rbin;
        S{i_T,i_r}.gr=S_cur.gr;
        S{i_T,i_r}.G6=reshape(S_cur.Gcorrs_mean(4,i_r,:),1,[]);
        rmax_vals(i_r)=S_cur.rmax_vals(i_r);
    end
end

%%
for i_r = 1:N_rmax
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]); 
    c_map=turbo(N_T+1); c_map=c_map(2:end,:);
    for i_T=1:N_T
        dispname=sprintf("$T = %.2f$",T_vals(i_T));
        plot(S{i_T,i_r}.Psibin,S{i_T,i_r}.Psi,...
            'DisplayName',dispname,...
            'Color',c_map(i_T,:));
        hold on;
    end
    hLegend=legend('Location','northeast','Interpreter','latex',...
        'NumColumns',2);
    
    xlim([0 1]);
    hXLabel = xlabel('$|\Psi|$','interpreter','latex');
    hYLabel = ylabel('$P(|\Psi|)$','interpreter','latex');
    if rmax_vals(i_r) > 0
        titlestr=sprintf("$r_{\\textrm{NB}} = %.2f$",rmax_vals(i_r));
    else
        titlestr="Delaunay Triangulation";
    end
    htitle=title(titlestr,'interpreter','latex');
    h_axis = gca;
    
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel, htitle], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation);
    
    
    figname=sprintf('%s/Psi_histogram/mxy_Psihist_sqrtN_%d_rmax_%.2f',figbase,sqrtN,rmax_vals(i_r));
    fprintf('Creating figure %s\n',figname)
        
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

%%
for i_r = 1:N_rmax
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]); 
    c_map=turbo(N_T+1); c_map=c_map(2:end,:);
    for i_T=1:N_T
        dispname=sprintf("$T = %.2f$",T_vals(i_T));
        loglog(S{i_T,i_r}.rbin/r_trig,abs(S{i_T,i_r}.G6),..../S{i_T}.gr,...
            'DisplayName',dispname,...
            'Color',c_map(i_T,:));
        hold on;
    end
    hLegend=legend('Location','northeast','Interpreter','latex',...
        'NumColumns',2);
    
    xlim([.6 6]);
    ylim([1e-4 5e0]);
    hXLabel = xlabel('$r/r_{\textrm{trig}}$','interpreter','latex');
    hYLabel = ylabel('$|G_6(r)|$','interpreter','latex');
    if rmax_vals(i_r) > 0
        titlestr=sprintf("$r_{\\textrm{NB}} = %.2f$",rmax_vals(i_r));
    else
        titlestr="Delaunay Triangulation";
    end
    htitle=title(titlestr,'interpreter','latex');
    h_axis = gca;
    
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel, htitle], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation);
    
    set(gca,'XTick',[.1:.1:1,1.5,2:1:10,15,20],'YTick',[.01 .1 1.0]);
    
    figname=sprintf('%s/G6/mxy_G6_sqrtN_%d_rmax_%.2f',figbase,sqrtN,rmax_vals(i_r));
    fprintf('Creating figure %s\n',figname)
        
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end