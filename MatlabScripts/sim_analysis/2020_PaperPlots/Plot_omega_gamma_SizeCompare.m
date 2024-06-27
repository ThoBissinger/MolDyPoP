%% 0 Initialization
clear
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/omega_gamma_plots',fig_base);

%% 0 Preparing data
model="mxy";modeltag="\textbf{MXY}";modeldir="mxy_3.00"; T = .17; T_dir = "T_.17"; simdir='210715_LinearTimeSampling';
% model="fmxy";modeltag="\textbf{FMXY}";modeldir="fmxy"; T = .17; T_dir = "T_.17"; simdir='211201_LongerTime';
% model="xy";modeltag="\textbf{XY}";modeldir="xy_s"; T = .89; T_dir = "T_.89"; simdir='211201_LongerTime';
sqrtN_vals = [16 32 64 128 256];
L_vals = [9.25 18.5 37 74 128];
N_N=numel(sqrtN_vals);

q_list=[1,3,6,9,12,16,19];

dt_sim = .01;
% 220201_ReducedSmapleStepDeltat
% 211201_LongerTime
% 210715_LinearTimeSampling
linestyle=':';
linewidth=1;
markersize=7;

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
t_cell=cell(1,N_N);
q_cell=cell(1,N_N);
cf_cell=cell(1,N_N);
ft_cell=cell(1,N_N);
om_vals_cell=cell(1,N_N);
gamma_cell=cell(1,N_N);
omega_1_cell=cell(1,N_N);
ft_peak_cell=cell(1,N_N);
om_peak_cell=cell(1,N_N);
chi_fit_cell=cell(1,N_N);
chi_vals_cell=cell(1,N_N);

for i_N=1:N_N
    sqrtN=sqrtN_vals(i_N);
    file=sprintf('/data/scc/thobi/%s/%s/sqrtN_%d/%s/samp_Dynamics.mat',...
        simdir,modeldir,sqrtN,T_dir);
    S=load(file,'gmperpmperp','averaging_times','qbin','chimperpq_av');
    t=S.averaging_times;
    q_vals=S.qbin;
%     q_cell{i_N}=q_vals(1:end-1);
    n_q=numel(q_vals);
    chi_vals_cell{i_N}=S.chimperpq_av/sqrtN^2;
    cf_full=real(S.gmperpmperp)/sqrtN^2;
    n_t=numel(t);
    t=S.averaging_times;
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
        t=S.averaging_times;
        cf = cf_full(q_indices);
        
        c = fit_DampedOscillator_RealSpace(t,cf/cf(1),n_period,weightexp,'omega_1');
    
        gamma_cur=c(1);
        om_cur=c(2);
        chi_cur=cf(1);
        gamma_vals(i_q)=gamma_cur;
        omega_1_vals(i_q)=om_cur;
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
    t_cell{i_N}=t_subcell{q_select};
    q_cell{i_N}=q_vals(q_select);
    cf_cell{i_N}=cf_subcell{q_select};
    ft_cell{i_N}=ft_subcell{q_select};
    om_vals_cell{i_N}=om_vals_subcell{q_select};
    gamma_cell{i_N}=gamma_vals(q_select);
    omega_1_cell{i_N}=omega_1_vals(q_select);
    ft_peak_cell{i_N}=ft_peak_vals(q_select);
    om_peak_cell{i_N}=om_peak_vals(q_select);
    chi_fit_cell{i_N}=chi_fit_vals(q_select);
end


























%% 0.75 Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall settings, annotatiosn etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.8*columnwidth_cm]);






%% 1 Subplot top: omega
ax_top=axes('InnerPosition',[.2 .5 .75 .4]);
pos_plottop=get(gca,'Position');
hold on;


% dispnames={"$2\gamma^{-1}(q)\chi_{m\perp}(q)$",... "$\frac{2}{\gamma(q)}\chi_{m\perp}(q)$"
%     "$S_{m\perp}^{\textrm{max}}(q)$"};
% y_data={2./gamma_vals.*chi_sim,ft_peak_vals};
c_map=linspecer(N_N);
for i_N = 3:N_N
    sqrtN=sqrtN_vals(i_N);
    dispname=sprintf('$N = (%d)^2$',sqrtN);
    q_vals=q_cell{i_N};
    y=omega_1_cell{i_N};
    h_plot{i_N} = plot(q_vals,y,...
        'DisplayName',dispname,...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'Marker', markers{i_N},...
        'MarkerSize', markersize,...
        'HandleVisibility', 'on', ...
        'Color',c_map(i_N,:));
end

c_guess=mean(y./q_vals);
qq=linspace(.07,.14);
plot(qq,2*c_guess*qq,'--',...
    'HandleVisibility','off',...'$\omega = cq',...
    'LineWidth',linewidth,...
    'Color','k');
text(1.2*qq(1),1.5*2*c_guess*qq(1),'$\omega = cq$',...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','right',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_labels);

xlim([.035 1]);
ylim([.007 Inf]);
% ylim([2e1 1e6]);


hXLabel = xlabel('$q$','interpreter','latex');
% hYLabel = ylabel('$S_{m\perp}(q,\omega) / \chi_{m\perp}(q)$','interpreter','latex');
% hYLabel = ylabel('$\gamma,S^{\textrm{max}}$','interpreter','latex');
hYLabel = ylabel('$\omega_1$','interpreter','latex');
hLegend = legend('Location', 'South','interpreter','latex',...
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
    'YScale','log', 'XScale','log',...
    'XTick',[.01 .02 .04 .08 .16 .32 .64 1.28],...
    'YTick',[.01 .02 .04 .08 .16 .32 .64],...
    'XAxisLocation','top',...
    'LineWidth', .5)
% ylim([1e1 Inf]);


pos = get(gca,'Position');
dim=[pos(1)+pos(3)-x_border-labelbox_size(1), ...
    pos(2)+y_border, ...
    labelbox_size(1),...
    labelbox_size(2)];
annotation('textbox',dim, 'String', 'a', 'FitBoxToText','off',...
    'interpreter','latex',...
    'Units','normalized',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');


annotation_str = {modeltag,sprintf('$T = %.2f$',T)};
pos = get(gca,'InnerPosition');
dim=[pos(1)+.04, ...
    pos(2)+pos(4)-.1, ...
    .125,.09];
h_annot=annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');












    
%% 2 Subplot bottom: omega fit vs FFT

% subplot(2,1,2,'replace');
% inset_pos=[pos(1)+.005,...
%     pos(2)+.4*pos(4),...
%     .6*pos(3),...
%     .56*pos(4)];
% ax_inset = axes('OuterPosition',inset_pos);
ax_bot=axes('InnerPosition',[.2 .1 .75 .4]);
pos_plotbottom=get(gca,'Position');



markers={'s', 'd', 'o', '^','v'};
c_map=linspecer(N_N);
for i_N = 3:N_N
    sqrtN=sqrtN_vals(i_N);
    dispname=sprintf('$N = (%d)^2$',sqrtN);
    q_vals=q_cell{i_N};
    y=gamma_cell{i_N};
    h_plot{i_N} = plot(q_vals,y,...
        'DisplayName',dispname,...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'Marker', markers{i_N},...
        'MarkerSize', markersize,...
        'HandleVisibility', 'on', ...
        'Color',c_map(i_N,:));
    hold on;
end

% c_guess=mean(y./q_vals);
qq=linspace(.07,.14);
fitob=fit(q_vals(:),y(:),'a*x^sigma'); % y_fit
a_fit=fitob.a;
sigma_fit=fitob.sigma;
y_fit=2*a_fit*qq.^sigma_fit;
plot(qq,y_fit,'--',...
    'HandleVisibility','off',...'$\omega = cq',...
    'LineWidth',linewidth,...
    'Color','k');
text(1.4*qq(1),2*y_fit(1),sprintf('$\\gamma = aq^{%.2f}$',sigma_fit),...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','right',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_labels);
xlim([.035 1]);
ylim([.0002 Inf]);


hXLabel = xlabel('$q$','interpreter','latex');
hYLabel = ylabel('$\gamma$','interpreter','latex');
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'XTick',[.01 .02 .04 .08 .16 .32 .64 1.28],...
    'YTick',[.0005 .001 .002 .004 .008 .016 .032 .064 .128 .256],...
    'YScale','log', 'XScale','log',...
    'LineWidth', .5)
% ylim([1e1 Inf]);


pos = get(gca,'Position');
dim=[pos(1)+pos(3)-x_border-labelbox_size(1), ...
    pos(2)+y_border, ...
    labelbox_size(1),...
    labelbox_size(2)];
annotation('textbox',dim, 'String', 'b', 'FitBoxToText','off',...
    'interpreter','latex',...
    'Units','normalized',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');



%% Saving
figname=sprintf('%s/%s_omegagamma_Compare_Size_T_%.3f',basedir,model,T);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end