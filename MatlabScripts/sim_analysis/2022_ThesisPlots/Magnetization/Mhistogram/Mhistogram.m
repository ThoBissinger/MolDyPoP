%% initialization
clear
pathbase="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/Magnetization/Mhistogram";
cd(pathbase)
run('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/initialization_script.m')
% pathbase=pwd;
model="mxy";
dataset_id="red";
basedir=sprintf('%s/dataset_%s',pathbase,dataset_id);


a=1.58;
K=2.16;
b=.934;
s=.373;
Pifunc=@(y) K*(exp(b*(y-s) - exp(b*(y-s)))).^a;

T=.11;
if (model == "mxy")
    sqrtN_vals=[16 32 64 128];
    L_vals=[9.25 18.5 37 74 148];
    T_str_vec=[".11" ".14" ".165" ".17" ".175" ".18" ".185"];
%     T_str_vec=[".17"];
    N_T = numel(T_str_vec);
%     data_root="/data/scc/thobi/211201_LongerTime/mxy_3.00/";
    if dataset_id == "lin"
        sqrtN_dirs=sprintfc("/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_%d",sqrtN_vals);
    elseif dataset_id == "red"
        sqrtN_dirs=sprintfc("/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00/sqrtN_%d",sqrtN_vals);
    elseif dataset_id == "long"
        sqrtN_dirs=sprintfc("/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_%d",sqrtN_vals);
    end
elseif (model == "fmxy")
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
elseif (model == "xy_s")
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
c_map=linspecer(N_N);
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
    hXLabel = xlabel('$(m - \langle m \rangle)/\sigma$','interpreter','latex');
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
    
    figname=sprintf('%s/%s_Mhistogram_T_%s',basedir,model,T_str);
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end


%%
sqrtN_vals=[16 32 64 128];
T_str=[".14",".17"];
for i_N = 1:4
    for i_T = 1:2
        S_lin{i_N,i_T}=load("/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_" + sqrtN_vals(i_N) + "/T_" + T_str(i_T) + "/run_1/output/sampling_output_Dynamics.mat",'absM','averaging_times'); 
        S_red{i_N,i_T}=load("/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00/sqrtN_" + sqrtN_vals(i_N) + "/T_" + T_str(i_T) + "/run_1/output/sampling_output_Dynamics.mat",'absM','averaging_times'); 
        S_long{i_N,i_T}=load("/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_" + sqrtN_vals(i_N) + "/T_" + T_str(i_T) + "/run_1/output/sampling_output_Dynamics.mat",'absM','averaging_times'); 
    
        n_samps(i_N,i_T,:) = [numel(S_lin{i_N,i_T}.absM),numel(S_red{i_N,i_T}.absM),numel(S_long{i_N,i_T}.absM)];
        tmax_vals(i_N,i_T,:) = [max(S_lin{i_N,i_T}.averaging_times),max(S_red{i_N,i_T}.averaging_times),max(S_long{i_N,i_T}.averaging_times)];
    end
end

%% Plot of figure parameters
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm],...
        'Resize','off');
loglog(sqrtN_vals, tmax_vals(:,1,3),'--^r',...
    'DisplayName',"Panel $(a)$");
hold on;
plot(sqrtN_vals, tmax_vals(:,2,3),'--vb',...
    'DisplayName',"Panel $(b)$");
plot(sqrtN_vals, tmax_vals(:,1,1),'--*g',...
    'DisplayName',"Panel $(c)$");
xlim([15,135])
ylim([400 30000])
h_axis = gca;
hXLabel = xlabel('$\sqrt{N}$','interpreter','latex');
hYLabel = ylabel('$t_{\textrm{max}}$','interpreter','latex');

set(gca,"XTick",sqrtN_vals,"XTickLabel",sqrtN_vals,...
    "YTick",[100 200 500 1000 2000 5000 10000 20000 50000]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legend, axes etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend = legend('Location', 'SouthEast','interpreter','latex','FontSize', fontsize_annotation,...
    'NumColumns',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Font
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca, 'FontName', 'cmr12','FontSize', fontsize_axis)
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels)

figname=sprintf('%s/%s_tmax_overview',pathbase,model);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');