%% 0 Initialization
clear
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/omega_gamma_plots',fig_base);

%% 0 Preparing data
model="mxy";modeltag="\textbf{MXY}";modeldir="mxy_3.00"; simdir='210715_LinearTimeSampling';
% model="fmxy";modeltag="\textbf{FMXY}";modeldir="fmxy"; T = .17; T_dir = "T_.17"; simdir='211201_LongerTime';
% model="xy";modeltag="\textbf{XY}";modeldir="xy_s"; T = .89; T_dir = "T_.89"; simdir='211201_LongerTime';
sqrtN_vals = [16 32 64 128 256];
L_vals = [9.25 18.5 37 74 148];
rho_vals = sqrtN_vals.^2 ./ L_vals.^2;
sqrtN_select=4:5;
N_N=numel(sqrtN_vals);


q_list=[1,3,6,9,12,16,19];
T_dirs = {"T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.25"};
T_vals = [.03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .25];
N_T = numel(T_vals);
simdirs_lowT={"220201_ReducedSmapleStepDeltat" "211201_LongerTime" "210715_LinearTimeSampling" ;...
    "220201_ReducedSmapleStepDeltat" "211201_LongerTime" "210715_LinearTimeSampling" ;...
    "220201_ReducedSmapleStepDeltat" "211201_LongerTime" "210715_LinearTimeSampling" ;...
    "220201_ReducedSmapleStepDeltat" "211201_LongerTime" "210715_LinearTimeSampling" ;...
    "220201_ReducedSmapleStepDeltat" "211201_LongerTime" "210715_LinearTimeSampling" };

simdirs_highT={"211201_LongerTime" "220201_ReducedSmapleStepDeltat" "210715_LinearTimeSampling";...
    "211201_LongerTime" "220201_ReducedSmapleStepDeltat" "210715_LinearTimeSampling";...
    "211201_LongerTime" "220201_ReducedSmapleStepDeltat" "210715_LinearTimeSampling" ;...
    "220201_ReducedSmapleStepDeltat" "211201_LongerTime" "210715_LinearTimeSampling" ;...
    "220201_ReducedSmapleStepDeltat" "211201_LongerTime" "210715_LinearTimeSampling" };

simdir_helicity = "210715_LinearTimeSampling";
T_sep=[.15 .15 .15 .15 .15];
dt_sim = .01;
% 220201_ReducedSmapleStepDeltat
% 211201_LongerTime
% 210715_LinearTimeSampling
linestyle='-';
linewidth=1;
markersize=5;

% colors={'r', 'b'};
markers={'s', 'd', 'o', '^','v'};


y_border = .008;
x_border = .008;
labelbox_size = [.05 .03];

runmax=500;

res_factor=3/4*log(runmax);
res_function = resolution_Gauss;
% res_function=@(t,tau) resolution_Laplace_pleateau(t,tau,tau);

%% 0.5 Data assembly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Peak Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_period = 4;   % number of periods included in DO fit
t_cell=cell(N_N,N_T);
q_cell=cell(N_N,N_T);
cf_cell=cell(N_N,N_T);
ft_cell=cell(N_N,N_T);
om_vals_cell=cell(N_N,N_T);
gamma_cell=cell(N_N,N_T);
omega_1_cell=cell(N_N,N_T);
ft_peak_cell=cell(N_N,N_T);
om_peak_cell=cell(N_N,N_T);
chi_fit_cell=cell(N_N,N_T);
chi_vals_cell=cell(N_N,N_T);
H_x_vals=zeros(N_N,N_T);
H_y_vals=zeros(N_N,N_T);
I_x_2_vals=zeros(N_N,N_T);
I_y_2_vals=zeros(N_N,N_T);
Helicity_vals=zeros(N_N,N_T);
for i_N=sqrtN_select
    sqrtN=sqrtN_vals(i_N);
    for i_T=1:N_T
        T_dir=T_dirs{i_T};
        T=T_vals(i_T);
        fprintf('sqrtN = %3d T = %.3f\n',sqrtN,T)
        
        curfile='';
        i=1;
        while (~isfile(curfile))
            if T < T_sep(i_N)
                cursimdir=simdirs_lowT{i_N,i};
            else
                cursimdir=simdirs_highT{i_N,i};
            end
            curfile=sprintf('/data/scc/thobi/%s/%s/sqrtN_%d/%s/samp_Dynamics.mat',...
                cursimdir,modeldir,sqrtN,T_dir);
            i=i+1;
        end
        S=load(curfile,'gmperpmperp','averaging_times','qbin','chimperpq_av');
        t=S.averaging_times + dt_sim;
        q_vals=S.qbin;
    %     q_cell{i_N}=q_vals(1:end-1);
        n_q=numel(q_vals);
        chi_vals_cell{i_N}=S.chimperpq_av/sqrtN^2;
        cf_full=real(S.gmperpmperp)/sqrtN^2;
        n_t=numel(t);
%         t=S.averaging_times;
        t_subcell=cell(1,n_q-1);
        cf_subcell=cell(1,n_q-1);
        ft_subcell=cell(1,n_q-1);
        om_vals_subcell=cell(1,n_q-1);
        gamma_vals=zeros(1,n_q-1);
        omega_1_vals=zeros(1,n_q-1);
        ft_peak_vals=zeros(1,n_q-1);
        om_peak_vals=zeros(1,n_q-1);
        chi_fit_vals=zeros(1,n_q-1);
    
        q_indices_base=(1:n_q:n_t*n_q)-1;
        q_select=intersect(q_list,1:n_q);
        for i_q = q_select
            q_indices = q_indices_base + i_q;
%             t=S.averaging_times;
            cf = cf_full(q_indices);
            
            c = fit_DampedOscillator_RealSpace(t,cf/cf(1),n_period,weightexp,'omega_1');
        
            gamma_cur=c(1);
            om_cur=c(2);
            chi_cur=cf(1);
            gamma_vals(i_q)=gamma_cur;
            omega_1_vals(i_q)=abs(om_cur);
            chi_fit_vals(i_q)=chi_cur;
            
        
            tau=res_factor/c(1);
    %         while tau > max(t)
    %             dt=(t(end)-t(1))/(numel(t)-1);
    %             [t,cf] = continue_cf(t,cf/cf(1),2,max(3*dt,om_cur/8/pi));
    %             cf = cf * chi_cur;
    %         end
            t_subcell{i_q}=t;
            cf_subcell{i_q}=cf;
            res_vals=res_function(t,tau);
        
            [ft_vals,om_vals]=FT_correlation(t, (cf .* res_vals), 0);
            ft_vals=real(ft_vals);
            ft_subcell{i_q}=ft_vals;
            om_vals_subcell{i_q}=om_vals;
            [ft_max,i_max]=max(ft_vals);
            om_peak_vals(i_q)=abs(om_vals(i_max));
            ft_peak_vals(i_q)=ft_max;
        end
        t_cell{i_N,i_T}=t_subcell{q_select};
        q_cell{i_N,i_T}=q_vals(q_select);
        cf_cell{i_N,i_T}=cf_subcell{q_select};
        ft_cell{i_N,i_T}=ft_subcell{q_select};
        om_vals_cell{i_N,i_T}=om_vals_subcell{q_select};
        gamma_cell{i_N,i_T}=gamma_vals(q_select);
        omega_1_cell{i_N,i_T}=omega_1_vals(q_select);
        ft_peak_cell{i_N,i_T}=ft_peak_vals(q_select);
        om_peak_cell{i_N,i_T}=om_peak_vals(q_select);
        chi_fit_cell{i_N,i_T}=chi_fit_vals(q_select);

        filename_Helicity=sprintf('/data/scc/thobi/%s/%s/sqrtN_%d/%s/samp_Dynamics_helicity.mat',...
                simdir_helicity,modeldir,sqrtN,T_dir);
        S_helicity=load(filename_Helicity,...
                "H_x","H_y","I_x_2","I_y_2");
        H_x_vals(i_N,i_T)=S_helicity.H_x;
        H_y_vals(i_N,i_T)=S_helicity.H_y;
        I_x_2_vals(i_N,i_T)=S_helicity.I_x_2; 
        I_y_2_vals(i_N,i_T)=S_helicity.I_y_2;
        Helicity_vals(i_N,i_T)=1/2/L_vals(i_N)^2 * ...
            ( abs(H_x_vals(i_N,i_T) + H_y_vals(i_N,i_T)) - ...
            1/T_vals(i_T)*(I_x_2_vals(i_N,i_T) + I_y_2_vals(i_N,i_T)));
    end
end




%% Evaluation of spins
c_sw=zeros(i_N,i_T);
c_sw_fit=zeros(i_N,i_T);
sigma_sw_fit=zeros(i_N,i_T);
fitob_sw=cell(i_N,i_T);
a_damp_fit=zeros(i_N,i_T);
sigma_damp_fit=zeros(i_N,i_T);
fitob_damp=cell(i_N,i_T);
for i_N=sqrtN_select
    sqrtN=sqrtN_vals(i_N);
    for i_T=1:N_T
        T_dir=T_dirs{i_T};
        T=T_vals(i_T);
        omega_1_vals=omega_1_cell{i_N,i_T};
        q_vals=q_cell{i_N,i_T};
        gamma_vals=gamma_cell{i_N,i_T};

        c_sw(i_N,i_T)=omega_1_vals(1)/q_vals(1);
%         c_sw(i_N,i_T)=mean(omega_1_vals./q_vals);
        
        fitob_sw{i_N,i_T}=fit(q_vals(:),omega_1_vals(:),'c*x^sigma',...
            'StartPoint',[c_sw(i_N,i_T),1]);
        c_sw_fit(i_N,i_T)=fitob_sw{i_N,i_T}.c;
        sigma_sw_fit(i_N,i_T)=fitob_sw{i_N,i_T}.sigma;

        init_val=gamma_vals(end)/q_vals(end)^2;
        fitob_damp{i_N,i_T}=fit(q_vals(:),gamma_vals(:),'a*x^sigma',...
            'StartPoint',[init_val,2]);
        a_damp_fit(i_N,i_T)=fitob_damp{i_N,i_T}.a;
        sigma_damp_fit(i_N,i_T)=fitob_damp{i_N,i_T}.sigma;
    end
end






















%% 0.75 Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall settings, annotatiosn etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.8*columnwidth_cm]);











    
%% 2 Subplot bottom: omega fit vs FFT

% subplot(2,1,2,'replace');
% inset_pos=[pos(1)+.005,...
%     pos(2)+.4*pos(4),...
%     .6*pos(3),...
%     .56*pos(4)];
% ax_inset = axes('OuterPosition',inset_pos);
% ax_bot=axes('InnerPosition',[.2 .1 .75 .4]);
% pos_plotbottom=get(gca,'Position');
ax_bot=subplot(2,1,2,'replace');
set(ax_bot,'InnerPosition',[.2 .1 .75 .4]);
pos_plottop=get(gca,'Position');
hold on;


c_map=linspecer(N_N);
for i_N = sqrtN_select
    sqrtN=sqrtN_vals(i_N);
    dispname=sprintf('$N = (%d)^2$',sqrtN);
    y=sigma_damp_fit(i_N,:);
    h_plot{i_N} = plot(T_vals,y,...
        'DisplayName',dispname,...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'Marker', markers{i_N},...
        'MarkerSize', markersize,...
        'HandleVisibility', 'on', ...
        'Color',c_map(i_N,:));
end
xlim([0 .25]);
ylim([0 3]);


hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$\sigma_{\gamma}$','interpreter','latex');
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
    'XGrid', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','lin',...
    'XAxisLocation','bottom',...
    'LineWidth', .5)
% ylim([1e1 Inf]);


% pos = get(gca,'Position');
% dim=[pos(1)+pos(3)-x_border-labelbox_size(1), ...
%     pos(2)+y_border, ...
%     labelbox_size(1),...
%     labelbox_size(2)];
% textpos=[.94 .06];
% text('Position',textpos, 'String', '(b)', ...
%     'interpreter','latex',...
%     'Units','normalized',...
%     'Color','black','FontSize', fontsize_subfiglabels);
h_text=add_subfig_label(gca,"$(b)$","nw","lin","lin",fontsize_subfiglabels);






%% 1 Subplot top: omega
ax_top=subplot(2,1,1,'replace');
set(ax_top,'InnerPosition',[.2 .5 .75 .4]);
pos_plottop=get(gca,'Position');
hold on;


c_map=linspecer(N_N);
for i_N = sqrtN_select
    sqrtN=sqrtN_vals(i_N);
    dispname=sprintf('$N = (%d)^2$',sqrtN);
%     y=c_sw_fit(i_N,:);
    y=c_sw(i_N,:);
    h_plot{i_N} = plot(T_vals,y,...
        'DisplayName',dispname,...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'Marker', markers{i_N},...
        'MarkerSize', markersize,...
        'HandleVisibility', 'on', ...
        'Color',c_map(i_N,:));
end
xlim([0 .25]);
ylim([0 .3]);


hXLabel = xlabel('$T$','interpreter','latex');
% hYLabel = ylabel('$S_{m\perp}(q,\omega) / \chi_{m\perp}(q)$','interpreter','latex');
% hYLabel = ylabel('$\gamma,S^{\textrm{max}}$','interpreter','latex');
hYLabel = ylabel('$c$','interpreter','latex');
% hLegend = legend('Location', 'SouthWest','interpreter','latex',...
%     'NumColumns',1);
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
% set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','lin',...
    'XAxisLocation','top',...
    'LineWidth', .5)
% ylim([1e1 Inf]);


% pos = get(gca,'Position');
% dim=[pos(1)+pos(3)-x_border-labelbox_size(1), ...
%     pos(2)+y_border, ...
%     labelbox_size(1),...
%     labelbox_size(2)];
% textpos=[.94 .06];
% text('Position',textpos, 'String', '(a)', ...
%     'interpreter','latex',...
%     'Units','normalized',...
%     'Color','black','FontSize', fontsize_subfiglabels);
h_text=add_subfig_label(gca,"$(a)$","nw","lin","lin",fontsize_subfiglabels);


annotation_str = {modeltag};
pos = get(gca,'InnerPosition');
dim=[pos(1)+pos(3)-.14, ...
    pos(2)+pos(4)-.075, ...
    .125,.09];
h_annot=annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=get(gca,"InnerPosition");
axes("Position",[pos(1)+.08, pos(2)+.03, pos(3)/2, pos(4)/2])

for i_N = sqrtN_select
    sqrtN=sqrtN_vals(i_N);
    dispname=sprintf('$N = (%d)^2$',sqrtN);
    h_plot{i_N} = plot(T_vals,rho_vals(i_N)*c_sw(i_N,:).^2./Helicity_vals(i_N,:),...
        'DisplayName',dispname,...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'Marker', 'o',...
        'MarkerSize', markersize,...
        'HandleVisibility', 'on', ...
        'Color',c_map(i_N,:));
    hold on;
end
xlims=[0 .19];
ylims=[.98 1.02];
xlim(xlims);
ylim(ylims);
text(xlims(1)+0.01,ylims(2)-.002,'$\frac{c^2}{\Upsilon / \rho}$',...
        'VerticalAlignment','top','HorizontalAlignment','left',...
        'interpreter','latex',...
        'Color','k','FontSize', fontsize_annotation);
text(xlims(2)-0.01,ylims(1)+.002,'$T$',...
        'VerticalAlignment','bottom','HorizontalAlignment','right',...
        'interpreter','latex',...
        'Color','k','FontSize', fontsize_annotation);    
set(gca, 'FontName', 'cmr12','FontSize', fontsize_axis,...
    'YTick',ylims(1):.01:ylims(2));







%% Saving
figname=sprintf('%s/%s_omegagamma_SpinwaveGammaexp',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end

%% 
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
for i_N = 1:N_N
    sqrtN=sqrtN_vals(i_N);
    dispname=sprintf('$N = (%d)^2$',sqrtN);
%     y=c_sw_fit(i_N,:);
    h_plot{i_N} = plot(T_vals,c_sw(i_N,:).^2./Helicity_vals(i_N,:),...
        'DisplayName',dispname,...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'Marker', 'o',...
        'MarkerSize', markersize,...
        'HandleVisibility', 'on', ...
        'Color',c_map(i_N,:));
    hold on;
%     h_plot{i_N} = plot(T_vals,Helicity_vals(i_N,:),...
%         'DisplayName',dispname,...
%         'LineStyle', linestyle, ...
%         'LineWidth', linewidth, ...
%         'Marker', 'x',...
%         'MarkerSize', markersize,...
%         'HandleVisibility', 'on', ...
%         'Color',c_map(i_N,:));
end
xlim([0 .18])
% ylim([0 2])
hXLabel = xlabel('$T$','interpreter','latex');
% hYLabel = ylabel('$S_{m\perp}(q,\omega) / \chi_{m\perp}(q)$','interpreter','latex');
% hYLabel = ylabel('$\gamma,S^{\textrm{max}}$','interpreter','latex');
hYLabel = ylabel('$c^2/\Upsilon$','interpreter','latex');
hLegend = legend('Location', 'NorthWest','interpreter','latex',...
    'NumColumns',1);
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation);

figname=sprintf('%s/%s_omegagamma_c_by_Upsilon',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end
