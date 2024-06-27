%% 0 Initialization mxy
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/Mhistogram',fig_base);
% fontsize_axis=15;
% fontsize_annotation=15;

% model_type = "mxy";
% model_type = "mxy";
% model_type = "fmxy";
model_type = "mxy";



a=1.58;
K=2.16;
b=.934;
s=.373;
Pifunc=@(y) K*(exp(b*(y-s) - exp(b*(y-s)))).^a;

T=.11;
if (model_type == "mxy")
    sqrtN_vals=[16 32 64 128];
    L_vals=[9.25 18.5 37 74 148];
    T_str_vec=[".11" ".14" ".165" ".17" ".175" ".18" ".185"];
%     T_str_vec=[".17"];
    N_T = numel(T_str_vec);
%     data_root="/data/scc/thobi/211201_LongerTime/mxy_3.00/";
    sqrtN_dirs=sprintfc("/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_%d",sqrtN_vals);
    sqrtN_dirs=sprintfc("/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00/sqrtN_%d",sqrtN_vals);
    sqrtN_dirs=sprintfc("/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_%d",sqrtN_vals);
%     sqrtN_dirs=["/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_16"
%         "/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_32"
%         "/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_64"
%         "/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_128"
%         "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256"
%         ];
elseif (model_type == "fmxy")
    sqrtN_vals=[16 32 64 128];
    L_vals=[9.25 18.5 37 74 148];
    T_str_vec=[".11" ".14" ".165" ".17" ".175" ".18" ".185"];
    N_T = numel(T_str_vec);
%     data_root="/data/scc/thobi/211201_LongerTime/mxy_3.00/";
    sqrtN_dirs=sprintfc("/data/scc/thobi/211201_LongerTime/fmxy/sqrtN_%d",sqrtN_vals);
%     sqrtN_dirs=["/data/scc/thobi/211201_LongerTime/fmxy/sqrtN_16"
%     "/data/scc/thobi/211201_LongerTime/fmxy/sqrtN_32"
%     "/data/scc/thobi/211201_LongerTime/fmxy/sqrtN_64"
%     "/data/scc/thobi/211201_LongerTime/fmxy/sqrtN_128"
%     "/data/scc/thobi/210715_LinearTimeSampling/fmxy/sqrtN_256"
%     ];
elseif (model_type == "xy_s")
    sqrtN_vals=[16 32 64 128];
    L_vals=sqrtN_vals;
    T_str_vec=[".85" ".91" ".95" "1.00"];
    N_T = numel(T_str_vec);
    sqrtN_dirs=sprintfc("/data/scc/thobi/211201_LongerTime/xy_s/sqrtN_%d",sqrtN_vals);
%     sqrtN_dirs=sprintfc("/data/scc/thobi/220201_ReducedSmapleStepDeltat/xy_s/sqrtN_%d",sqrtN_vals);
end
symbols_vec=['o','v','^','s','d'];
sampfilename="samp_Dynamics_Mhistogram";
N_N=numel(sqrtN_vals);
%% 1 Plotting the data collapse
c_map=colormap_sqrtN();
for i_T = 1:N_T
    T_str = T_str_vec(i_T);
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm],...
        'Resize','off');
    for i_N = 1:N_N
        sqrtN = sqrtN_vals(i_N);
        dispname = sprintf('$N = (%d)^2$',sqrtN);
        sampfile_cur=sprintf('%s/T_%s/%s',sqrtN_dirs{i_N},T_str,sampfilename);
        load(sampfile_cur);
%         semilogy(M_edges,sqrt(M_var)*M_pdf,...
        semilogy(M_edges,M_pdf,...
            'DisplayName',dispname, ...
            'Color',c_map(i_N,:),...
            'LineStyle','none',...
            'LineWidth',1,...
            'Marker',symbols_vec(i_N),...
            'MarkerSize',3,...
            'MarkerIndices',[1:2:numel(M_edges)]);
        hold on;
    end
    semilogy(M_edges,Pifunc(M_edges),'--',...
        'Color','black',...
        'LineWidth',2,...
        'DisplayName','universal $\Pi$');
    xlim([-6,3])
    h_axis = gca;
    hXLabel = xlabel('$(m - \langle m \rangle)/\sigma_m$','interpreter','latex');
    hYLabel = ylabel('$P(m)$','interpreter','latex');
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Legend, axes etc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hLegend = legend('Location', 'South','interpreter','latex','FontSize', fontsize_annotation,...
        'NumColumns',2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Font
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(gca, 'FontName', 'cmr12','FontSize', fontsize_axis)
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels)
    
    figname=sprintf('%s/%s_Mhistogram_T_%s',basedir,model_type,T_str);
    %                 figname=sprintf('%s/SbyL_vsomL/%s_S_%s_FS_n_q_%d_z_%.2f_tau_%de3_T_%.3f',basedir,curmodel,plot_type_cur,n_q,z,tau_laplace/1e3,T);
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end