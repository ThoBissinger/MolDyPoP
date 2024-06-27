%% 0 Initialization
clear
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/TCF/Multiplot',fig_base);

%% 0 Preparing data
sqrtN_vals = [16 32 64 128];
L_vals = [9.25 18.5 37 74];

sqrtN=128;

T_vals = [.11 .14 .17 .185];
T_dirs = {"T_.11" "T_.14" "T_.17" "T_.185"};
N_T = numel(T_vals);

i_q = 3;

mxy_dir=sprintf("/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_%d",sqrtN);

fmxy_dir=sprintf("/data/scc/thobi/211201_LongerTime/fmxy/sqrtN_%d",sqrtN);
% fmxy_dir_shorttime="/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00";

mxy_runmax=500;
fmxy_runmax=500;

mxy_res_factor=3/4*log(mxy_runmax);
fmxy_res_factor=3/4*log(fmxy_runmax);
res_function = resolution_Gauss;

marker_types=["v" "^" "o" "d"];

mxy_mperp_cf=cell(size(T_vals));
mxy_qbin=cell(size(T_vals));
mxy_t=cell(size(T_vals));
fmxy_mperp_cf=cell(size(T_vals));
fmxy_qbin=cell(size(T_vals));
fmxy_t=cell(size(T_vals));

for i_T = 1:N_T
    mxy_curfile=sprintf('%s/%s/samp_Dynamics_collect.mat',mxy_dir,T_dirs{i_T});
    fmxy_curfile=sprintf('%s/%s/samp_Dynamics_collect.mat',fmxy_dir,T_dirs{i_T});

    load(mxy_curfile,'averaging_times','qbin','gmperpmperp_collect',...
        'gww_collect','gxx_collect','gyy_collect');
    mxy_mperp_cf{i_T} = real(gmperpmperp_collect);
    mxy_w_cf{i_T} = real(gww_collect);    
    mxy_m_cf{i_T} = real(gxx_collect) + real(gyy_collect);        
    mxy_qbin{i_T} = qbin;
    mxy_t{i_T} = averaging_times;

    load(fmxy_curfile,'averaging_times','qbin','gmperpmperp_collect',...
        'gww_collect','gxx_collect','gyy_collect');
    fmxy_mperp_cf{i_T} = real(gmperpmperp_collect);
    fmxy_w_cf{i_T} = real(gww_collect);    
    fmxy_m_cf{i_T} = real(gxx_collect) + real(gyy_collect);        
    fmxy_qbin{i_T} = qbin;
    fmxy_t{i_T} = averaging_times;
end

%% 1 Subplot 1: mperp cf vs t for different T, MXY MODEL
subplot(3,2,1,'replace');

mxy_gamma_vals=zeros(size(T_vals));
mxy_omega_1_vals=zeros(size(T_vals));
mxy_chi_vals=zeros(size(T_vals));
zeros(size(T_vals));
n_period = 4;
t_max = 2e3;

c_map=linspecer(N_T);
for i_T = 1 : N_T
    T=T_vals(i_T);
    t=mxy_t{i_T};
    n_t=numel(t);
    qbin_vals=mxy_qbin{i_T};
    n_q = numel(qbin_vals);
    q = qbin_vals(i_q);
    
    dispname=sprintf('$T = %.3f$', T);
    q_indices = (i_q):n_q:n_q*n_t;
    cf = mxy_mperp_cf{i_T}(:,q_indices);
    cf_mean = mean(cf);
    mxy_chi_vals(i_T) = cf_mean(1);
    cf_mean = cf_mean / mxy_chi_vals(i_T);
    cf_std = std(cf / mxy_chi_vals(i_T));
    cf_err = cf_std/sqrt(mxy_runmax);
    
    c = fit_DampedOscillator_RealSpace(t,cf_mean/cf_mean(1),n_period,weightexp,'omega_1');
    mxy_gamma_vals(i_T) = c(1);
    mxy_omega_1_vals(i_T) = c(2);
    
%     cf= real(cf)/real(cf(1));
    offset= -2.5*(i_T - 1);
    
    yline(offset,'--',...
        'HandleVisibility','off');
    hold on;
    
    h_fit_data(i_T) = plot(t,fitfunc_DO(t,1,c) + offset,...
        'LineStyle', '-', ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:));
    
%     err_indices=1:4:numel(cf_mean);
%     h_sim_data(i_T) = errorbar(t(err_indices),cf_mean(err_indices) + offset, cf_err(err_indices),...
%         'LineStyle', 'none', ...
%         'HandleVisibility', 'off', ...
%         'Color',c_map(i_T,:));
    h_sim_data(i_T) = plot(t,cf_mean + offset,...
        'LineStyle', 'none', ...
        'Marker', marker_types(i_T), 'MarkerSize', 4, ...
        'LineWidth', .5, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:),...
        'MarkerFaceColor','none','MarkerEdgeColor',c_map(i_T,:), ...
        'MarkerIndices',1:15:numel(t));
%         'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k', ...

    h_leg_data(i_T) = plot(NaN*t,NaN*t,...
        'LineStyle', '-', ...
        'LineWidth', 1, ...
        'Marker', marker_types(i_T), 'MarkerSize', 4, ...
        'DisplayName', dispname, ...
        'MarkerFaceColor','none','MarkerEdgeColor',c_map(i_T,:), ...
        'Color',c_map(i_T,:));
    
end
ylim([offset-1.2 5])
xlim([0 t_max]);
hLegend = legend('Location', 'NorthEast','interpreter','latex',...
    'NumColumns',1);

hXLabel = xlabel('$t$','interpreter','latex');
hYLabel = ylabel('$C_{m\perp}(q,t) / \chi_{m\perp}(q)$','interpreter','latex');
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', [], ...
    'LineWidth', .5)


annotation_str = {sprintf('$N = (%d)^2$',sqrtN),sprintf('$q = %.3f$',q)};
dim=[.14 .872 .2 .05];
annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

dim=[.427 .703 .2 .05];
annotation('textbox',dim, 'String', 'a', 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');



















%% 2 Subplot 2: mperp cf vs t for different T, FMXY MODEL
subplot(3,2,2,'replace');

fmxy_gamma_vals=zeros(size(T_vals));
fmxy_omega_1_vals=zeros(size(T_vals));
fmxy_chi_vals=zeros(size(T_vals));
zeros(size(T_vals));
n_period = 4;
t_max = 2e3;

c_map=linspecer(N_T);
for i_T = 1 : N_T
    T=T_vals(i_T);
    t=fmxy_t{i_T};
    n_t=numel(t);
    qbin_vals=fmxy_qbin{i_T};
    n_q = numel(qbin_vals);
    q = qbin_vals(i_q);
    
    dispname=sprintf('$T = %.3f$', T);
    q_indices = (i_q):n_q:n_q*n_t;
    cf = fmxy_mperp_cf{i_T}(:,q_indices);
    cf_mean = mean(cf);
    fmxy_chi_vals(i_T) = cf_mean(1);
    cf_mean = cf_mean / fmxy_chi_vals(i_T);
    cf_std = std(cf / fmxy_chi_vals(i_T));
    cf_err = cf_std/sqrt(fmxy_runmax);
    
    c = fit_DampedOscillator_RealSpace(t,cf_mean/cf_mean(1),n_period,weightexp,'omega_1');
    fmxy_gamma_vals(i_T) = c(1);
    fmxy_omega_1_vals(i_T) = c(2);
    
    offset= -2.5*(i_T - 1);
    
    yline(offset,'--',...
        'HandleVisibility','off');
    hold on;
    
    h_fit_data(i_T) = plot(t,fitfunc_DO(t,1,c) + offset,...
        'LineStyle', '-', ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:));
    
    h_sim_data(i_T) = plot(t,cf_mean + offset,...
        'LineStyle', 'none', ...
        'Marker', marker_types(i_T), 'MarkerSize', 4, ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:),...
        'MarkerFaceColor','none','MarkerEdgeColor',c_map(i_T,:), ...
        'MarkerIndices',1:4:numel(t));
    
end
ylim([offset-1.2 5])
xlim([0 t_max]);

hXLabel = xlabel('$t$','interpreter','latex');
hYLabel = ylabel('$C_{m\perp}(q,t) / \chi_{m\perp}(q)$','interpreter','latex');
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', [], ...
    'LineWidth', .5)



dim=[.868 .703 .2 .05];
annotation('textbox',dim, 'String', 'b', 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');





















%% 3 Subplot 3: mperp cf FFT vs omega for different T, MXY MODEL
subplot(3,2,3,'replace');

mxy_gamma_vals=zeros(size(T_vals));
mxy_omega_1_vals=zeros(size(T_vals));
mxy_chi_vals=zeros(size(T_vals));
om_max_vals_mxy=zeros(size(T_vals));
ft_max_vals_mxy=zeros(size(T_vals));
ft_zero_vals_mxy=zeros(size(T_vals));
ft_sidepeak_max_mxy=zeros(size(T_vals));

n_period = 4;
t_max = 2e3;

c_map=linspecer(N_T);
for i_T = 1 : N_T
    T=T_vals(i_T);
    t=mxy_t{i_T};
    n_t=numel(t);
    qbin_vals=mxy_qbin{i_T};
    n_q = numel(qbin_vals);
    q = qbin_vals(i_q);
    
    dispname=sprintf('$T = %.3f$', T);
    q_indices = (i_q):n_q:n_q*n_t;
    cf = mxy_mperp_cf{i_T}(:,q_indices);
    cf_mean = mean(cf);
    mxy_chi_vals(i_T) = cf_mean(1);
    c = fit_DampedOscillator_RealSpace(t,cf_mean/cf_mean(1),n_period,weightexp,'omega_1');
    mxy_gamma_vals(i_T) = c(1);
    mxy_omega_1_vals(i_T) = c(2);
    

    tau=mxy_res_factor/mxy_gamma_vals(i_T);
    res_vals=res_function(t,tau);

    [ft_vals,om_vals_cur]=FT_correlation(t, (cf .* res_vals)' /cf_mean(1), 0);
    ft_vals=real(ft_vals)';
    ft_mean=mean(real(ft_vals));
    ft_std=std(real(ft_vals));
    ft_err=ft_std/sqrt(mxy_runmax);
    ft_var=var(real(ft_vals));
    om_vals=om_vals_cur;

    [ft_max,i_max]=max(ft_mean);
    om_max_vals_mxy(i_T) = abs(om_vals_cur(i_max));
    ft_max_vals_mxy(i_T) = ft_max;

    i_zero = find(om_vals_cur >= 0, 1);
    ft_zero_vals_mxy(i_T) = ft_mean(i_zero);

    i_sidepeak = find(abs(om_vals_cur) < (om_max_vals_mxy(i_T) + mxy_gamma_vals(i_T))/2);
    ft_sidepeak_max_mxy(i_T) = max(ft_mean(i_sidepeak));
    
    
%     h_sim_data_line(i_T) = plot(om_vals,ft_mean,...
%         'LineStyle', '-', ...
%         'LineWidth', 1, ...
%         'HandleVisibility', 'off', ...
%         'Color','black');
%     hold on;
%     h_sim_data(i_T) = errorbar(om_vals(1:2:end),ft_mean(1:2:end),ft_err(1:2:end),...
%         'LineStyle', 'none', ...
%         'LineWidth', 1, ...
%         'HandleVisibility', 'off', ...
%         'Color',c_map(i_T,:));
    hold on;
    h_sim_data_line(i_T) = plot(om_vals,ft_mean,...
        'LineStyle', '-', ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color','black');
    h_sim_data(i_T) = errorbar(om_vals(1:3:end),ft_mean(1:3:end),ft_err(1:3:end),...
        'LineStyle', 'none', ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:));

end
ylim([0 Inf])
xlim([0 .06]);

hXLabel = xlabel('$\omega$','interpreter','latex');
hYLabel = ylabel('$S_{m\perp}(q,\omega) / \chi_{m\perp}(q)$','interpreter','latex');
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin',...
    'LineWidth', .5)



dim=[.426 .404 .2 .05];
annotation('textbox',dim, 'String', 'c', 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');










%% 4 Subplot 4: mperp cf FFT vs omega for different T, FMXY MODEL
subplot(3,2,4,'replace');

fmxy_gamma_vals=zeros(size(T_vals));
fmxy_omega_1_vals=zeros(size(T_vals));
fmxy_chi_vals=zeros(size(T_vals));
om_max_vals_fmxy=zeros(size(T_vals));
ft_max_vals_fmxy=zeros(size(T_vals));
ft_zero_vals_fmxy=zeros(size(T_vals));
ft_sidepeak_max_fmxy=zeros(size(T_vals));

n_period = 4;
t_max = 2e3;

c_map=linspecer(N_T);
for i_T = 1 : N_T
    T=T_vals(i_T);
    t=fmxy_t{i_T};
    n_t=numel(t);
    qbin_vals=fmxy_qbin{i_T};
    n_q = numel(qbin_vals);
    q = qbin_vals(i_q);
    
    dispname=sprintf('$T = %.3f$', T);
    q_indices = (i_q):n_q:n_q*n_t;
    cf = fmxy_mperp_cf{i_T}(:,q_indices);
    cf_mean = mean(cf);
    fmxy_chi_vals(i_T) = cf_mean(1);
    c = fit_DampedOscillator_RealSpace(t,cf_mean/cf_mean(1),n_period,weightexp,'omega_1');
    fmxy_gamma_vals(i_T) = c(1);
    fmxy_omega_1_vals(i_T) = c(2);
    

    tau=fmxy_res_factor/fmxy_gamma_vals(i_T);
    res_vals=res_function(t,tau);

    [ft_vals,om_vals_cur]=FT_correlation(t, (cf .* res_vals)' /cf_mean(1), 0);
    ft_vals=real(ft_vals)';
    ft_mean=mean(real(ft_vals));
    ft_std=std(real(ft_vals));
    ft_err=ft_std/sqrt(fmxy_runmax);
    ft_var=var(real(ft_vals));
    om_vals=om_vals_cur;

    [ft_max,i_max]=max(ft_mean);
    om_max_vals_fmxy(i_T) = abs(om_vals_cur(i_max));
    ft_max_vals_fmxy(i_T) = ft_max;

    i_zero = find(om_vals_cur >= 0, 1);
    ft_zero_vals_fmxy(i_T) = ft_mean(i_zero);

    i_sidepeak = find(abs(om_vals_cur) < (om_max_vals_fmxy(i_T) + fmxy_gamma_vals(i_T))/2);
    ft_sidepeak_max_fmxy(i_T) = max(ft_mean(i_sidepeak));
    
    
%     h_sim_data(i_T) = errorbar(om_vals,ft_mean,ft_err,...
%         'LineStyle', '-', ...
%         'HandleVisibility', 'off', ...
%         'Color',c_map(i_T,:));
    hold on;
    h_sim_data_line(i_T) = plot(om_vals,ft_mean,...
        'LineStyle', '-', ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color','black');
    h_sim_data(i_T) = errorbar(om_vals(1:1:end),ft_mean(1:1:end),ft_err(1:1:end),...
        'LineStyle', 'none', ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:));
end
ylim([0 Inf])
xlim([0 .06]);

hXLabel = xlabel('$\omega$','interpreter','latex');
hYLabel = ylabel('$S_{m\perp}(q,\omega) / \chi_{m\perp}(q)$','interpreter','latex');
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin',...
    'LineWidth', .5)



dim=[.868 .403 .2 .05];
annotation('textbox',dim, 'String', 'd', 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');























%% 5 Subplot 5: m cf FFT vs omega for different T, MXY MODEL
subplot(3,2,5,'replace');

mxy_gamma_vals=zeros(size(T_vals));
mxy_omega_1_vals=zeros(size(T_vals));
mxy_chi_vals=zeros(size(T_vals));
om_max_vals_mxy=zeros(size(T_vals));
ft_max_vals_mxy=zeros(size(T_vals));
ft_zero_vals_mxy=zeros(size(T_vals));
ft_sidepeak_max_mxy=zeros(size(T_vals));

n_period = 4;
t_max = 2e3;

c_map=linspecer(N_T);
for i_T = 1 : N_T
    T=T_vals(i_T);
    t=mxy_t{i_T};
    n_t=numel(t);
    qbin_vals=mxy_qbin{i_T};
    n_q = numel(qbin_vals);
    q = qbin_vals(i_q);
    
    dispname=sprintf('$T = %.3f$', T);
    q_indices = (i_q):n_q:n_q*n_t;
    cf = mxy_m_cf{i_T}(:,q_indices);
    cf_mean = mean(cf);
    mxy_chi_vals(i_T) = cf_mean(1);
    c = fit_DampedOscillator_RealSpace(t,cf_mean/cf_mean(1),n_period,weightexp,'omega_1');
    mxy_gamma_vals(i_T) = c(1);
    mxy_omega_1_vals(i_T) = c(2);
    

    tau=mxy_res_factor/mxy_gamma_vals(i_T);
    res_vals=res_function(t,tau);

    [ft_vals,om_vals_cur]=FT_correlation(t, (cf .* res_vals)' /cf_mean(1), 0);
    ft_vals=real(ft_vals)';
    ft_mean=mean(real(ft_vals));
    ft_std=std(real(ft_vals));
    ft_err=ft_std/sqrt(mxy_runmax);
    ft_var=var(real(ft_vals));
    om_vals=om_vals_cur;

    [ft_max,i_max]=max(ft_mean);
    om_max_vals_mxy(i_T) = abs(om_vals_cur(i_max));
    ft_max_vals_mxy(i_T) = ft_max;

    i_zero = find(om_vals_cur >= 0, 1);
    ft_zero_vals_mxy(i_T) = ft_mean(i_zero);

    i_sidepeak = find(abs(om_vals_cur) < (om_max_vals_mxy(i_T) + mxy_gamma_vals(i_T))/2);
    ft_sidepeak_max_mxy(i_T) = max(ft_mean(i_sidepeak));
    
    
%     h_sim_data(i_T) = errorbar(om_vals,ft_mean,ft_err,...
%         'LineStyle', '-', ...
%         'HandleVisibility', 'off', ...
%         'Color',c_map(i_T,:));
    hold on;
    h_sim_data_line(i_T) = plot(om_vals,ft_mean,...
        'LineStyle', '-', ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color','black');
    h_sim_data(i_T) = errorbar(om_vals(1:3:end),ft_mean(1:3:end),ft_err(1:3:end),...
        'LineStyle', 'none', ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:));
    
end
ylim([0 Inf])
xlim([0 .06]);

hXLabel = xlabel('$\omega$','interpreter','latex');
hYLabel = ylabel('$S_{m}(q,\omega) / \chi_{m}(q)$','interpreter','latex');
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin',...
    'LineWidth', .5)



dim=[.426 .105 .2 .05];
annotation('textbox',dim, 'String', 'e', 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');
















%% 6 Subplot 6: m cf FFT vs omega for different T, FMXY MODEL
subplot(3,2,6,'replace');

fmxy_gamma_vals=zeros(size(T_vals));
fmxy_omega_1_vals=zeros(size(T_vals));
fmxy_chi_vals=zeros(size(T_vals));
om_max_vals_fmxy=zeros(size(T_vals));
ft_max_vals_fmxy=zeros(size(T_vals));
ft_zero_vals_fmxy=zeros(size(T_vals));
ft_sidepeak_max_fmxy=zeros(size(T_vals));

n_period = 4;
t_max = 2e3;

c_map=linspecer(N_T);
for i_T = 1 : N_T
    T=T_vals(i_T);
    t=fmxy_t{i_T};
    n_t=numel(t);
    qbin_vals=fmxy_qbin{i_T};
    n_q = numel(qbin_vals);
    q = qbin_vals(i_q);
    
    dispname=sprintf('$T = %.3f$', T);
    q_indices = (i_q):n_q:n_q*n_t;
    cf = fmxy_m_cf{i_T}(:,q_indices);
    cf_mean = mean(cf);
    fmxy_chi_vals(i_T) = cf_mean(1);
    c = fit_DampedOscillator_RealSpace(t,cf_mean/cf_mean(1),n_period,weightexp,'omega_1');
    fmxy_gamma_vals(i_T) = c(1);
    fmxy_omega_1_vals(i_T) = c(2);
    

    tau=fmxy_res_factor/fmxy_gamma_vals(i_T);
    res_vals=res_function(t,tau);

    [ft_vals,om_vals_cur]=FT_correlation(t, (cf .* res_vals)' /cf_mean(1), 0);
    ft_vals=real(ft_vals)';
    ft_mean=mean(real(ft_vals));
    ft_std=std(real(ft_vals));
    ft_err=ft_std/sqrt(fmxy_runmax);
    ft_var=var(real(ft_vals));
    om_vals=om_vals_cur;

    [ft_max,i_max]=max(ft_mean);
    om_max_vals_fmxy(i_T) = abs(om_vals_cur(i_max));
    ft_max_vals_fmxy(i_T) = ft_max;

    i_zero = find(om_vals_cur >= 0, 1);
    ft_zero_vals_fmxy(i_T) = ft_mean(i_zero);

    i_sidepeak = find(abs(om_vals_cur) < (om_max_vals_fmxy(i_T) + fmxy_gamma_vals(i_T))/2);
    ft_sidepeak_max_fmxy(i_T) = max(ft_mean(i_sidepeak));
    
    
%     h_sim_data(i_T) = errorbar(om_vals,ft_mean,ft_err,...
%         'LineStyle', '-', ...
%         'HandleVisibility', 'off', ...
%         'Color',c_map(i_T,:));
    hold on;
    h_sim_data_line(i_T) = plot(om_vals,ft_mean,...
        'LineStyle', '-', ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color','black');
    h_sim_data(i_T) = errorbar(om_vals(1:1:end),ft_mean(1:1:end),ft_err(1:1:end),...
        'LineStyle', 'none', ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:));
    
end
ylim([0 Inf])
xlim([0 .06]);

hXLabel = xlabel('$\omega$','interpreter','latex');
hYLabel = ylabel('$S_{m}(q,\omega) / \chi_{m}(q)$','interpreter','latex');
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin',...
    'LineWidth', .5)



dim=[.868 .105 .2 .05];
annotation('textbox',dim, 'String', 'f', 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall settings, annotatiosn etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_mxy_title=annotation('textbox',[.25 .92 .2 .05],...
    'string','MXY model','Units','normalized',...
    'FitBoxToText','on',...
    'EdgeColor','none',...
    'FontSize',fontsize_labels);

h_fmxy_title=annotation('textbox',[.67 .92 .2 .05],...
    'string','FMXY model','Units','normalized',...
    'FitBoxToText','on',...
    'EdgeColor','none',...
    'FontSize',fontsize_labels);

% set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
set(gcf,'units','centimeters','OuterPosition',[0 0 pagewidth_cm pageheight_cm]);
% set_fonts_default;























%% Saving
figname=sprintf('%s/Compare_Corrs_sqrtN_%d_q_%.3f',basedir,sqrtN,q);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end