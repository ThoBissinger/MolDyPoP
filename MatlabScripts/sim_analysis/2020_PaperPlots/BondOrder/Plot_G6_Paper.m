%% Initialization
clear all
cd /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/BondOrder
% addpath "/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/";
run ../initialization_script;
pathbase="/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00";
figbase="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/BondOrder";
cd(figbase);
sqrtN=256;
L=9.25*sqrtN/16;
rho=(sqrtN/L)^2;
T_str=[".01" ".03" ".05" ".09" ".14" ".20"];
% T_str=[ ".01"  ".03"  ".09"  ".15"  ".19"  ".31"];
T_vals=str2double(T_str);
N_T=numel(T_vals);
rmax_vals=[0 0.66 0.73 0.80 1.08];
N_rmax=numel(rmax_vals);
pathbase=pathbase + "/sqrtN_" + sqrtN;
i_r=2;
r_hex=sqrt(2/sqrt(3)/rho);
% r_hex=sqrt(1/rho);
r_hon=sqrt(2)*r_hex;

S=cell(N_T,1);
for i_T = 1:N_T
    curfile=pathbase+"/T_" + T_str(i_T) + "/samp_Dynamics_bondorder.mat";
    S_cur=load(curfile,"rbin","Gcorrs_mean","rmax_vals","gr");
        S{i_T}.rbin=S_cur.rbin;
        S{i_T}.gr=S_cur.gr;
        S{i_T}.G6=reshape(S_cur.Gcorrs_mean(4,:,:),[],numel(S_cur.rbin));
end
rmax_vals=S_cur.rmax_vals;

%%
i_r=3;
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]); 
c_map=turbo(N_T+1); c_map=c_map(2:end,:);
i_T=1;
dispname=sprintf("$T = %.2f$",T_vals(i_T));
loglog(S{i_T}.rbin/r_hex,abs(real(S{i_T}.G6(2,:))),..../S{i_T}.gr,...
    'DisplayName',dispname,...
    'LineWidth',1.3,...
    'Color',c_map(i_T,:));
hold on;
for i_T=2:N_T
    dispname=sprintf("$T = %.2f$",T_vals(i_T));
    loglog(S{i_T}.rbin/r_hex,abs(real(S{i_T}.G6(i_r,:))),..../S{i_T}.gr,...
        'DisplayName',dispname,...
        'LineWidth',1.3,...
        'Color',c_map(i_T,:));
end
hLegend=legend('Location','northeast','Interpreter','latex',...
    'NumColumns',2);

xlim([.6 4]);
ylim([1e-4 1e1]);
hXLabel = xlabel('$r/r_{\textrm{hex}}$','interpreter','latex');
hYLabel = ylabel('$|G_6(r)|$','interpreter','latex');
% if rmax_vals(i_r) > 0
%     titlestr=sprintf("$r_{\\textrm{NB}} = %.2f$",rmax_vals(i_r));
% else
%     titlestr="Delaunay Triangulation";
% end
% htitle=title(titlestr,'interpreter','latex');
h_axis = gca;

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation);

set(gca,'XTick',[.1:.1:1,1.5,2:1:10,15,20],'YTick',[.001 .01 .1 1.0 10]);

figname=sprintf('%s/G6/mxy_G6_sqrtN_%d',figbase,sqrtN);
fprintf('Creating figure %s\n',figname)
    
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
