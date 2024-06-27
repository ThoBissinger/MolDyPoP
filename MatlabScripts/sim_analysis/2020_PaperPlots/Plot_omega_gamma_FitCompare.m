%% 0 Initialization
clear
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/omega_gamma_plots',fig_base);

%% 0 Preparing data
model="mxy";modeltag="\textbf{MXY}";modeldir="mxy_3.00"; T = .17; T_dir = "T_.17"; simdir='210715_LinearTimeSampling'; % simdir='220201_ReducedSmapleStepDeltat';
% model="fmxy";modeltag="\textbf{FMXY}";modeldir="fmxy"; T = .17; T_dir = "T_.17"; simdir='210715_LinearTimeSampling';
% model="xy";modeltag="\textbf{XY}";modeldir="xy_s"; T = .89; T_dir = "T_.89"; simdir='211201_LongerTime';
sqrtN_vals = [16 32 64 128];
L_vals = [9.25 18.5 37 74];

q_switch=0.5;

runmax=250;
runsep=50;
n_runbin = 5;

sqrtN=128;
last_q_switch = 1; % If 0, all q values can be used, skips the last entry if 1
dt_sim = .01;
% 211201_LongerTime
% 210715_LinearTimeSampling
file=sprintf('/data/scc/thobi/%s/%s/sqrtN_%d/%s/samp_Dynamics.mat',...
    simdir,modeldir,sqrtN,T_dir);
fullfile=sprintf('/data/scc/thobi/%s/%s/sqrtN_%d/%s/samp_Dynamics_collect.mat',...
    simdir,modeldir,sqrtN,T_dir);
S=load(file,'gmperpmperp','averaging_times','qbin','chimperpq_av');
load(fullfile,'gmperpmperp_collect');

% S=load(file);
t=S.averaging_times;
% q_vals=S.qbin(1:end-1);
q_vals=S.qbin;
chi_sim=S.chimperpq_av/sqrtN^2;
cf_full=real(S.gmperpmperp)/sqrtN^2;
chi_init=cf_full(1:numel(q_vals));
n_t=numel(t);
n_q=numel(q_vals);
q_vals = q_vals(1:end-last_q_switch);
chi_init = chi_init(1:end-last_q_switch);

y_border = .008;
x_border = .008;
labelbox_size = [.05 .03];

fullfile_betterRes=sprintf('/data/scc/thobi/220201_ReducedSmapleStepDeltat/%s/sqrtN_%d/%s/samp_Dynamics_collect.mat',...
    modeldir,sqrtN,T_dir);
if isfile(fullfile_betterRes)
    S_betterres=load(fullfile_betterRes,'gmperpmperp','averaging_times','qbin','gmperpmperp_collect');
    cf_full_betterres=real(S_betterres.gmperpmperp)/sqrtN^2;
    cf_collect_betterres=real(S_betterres.gmperpmperp_collect)/sqrtN^2;
    q_betterres=S_betterres.qbin;
end

res_factor=3/4*log(runmax);
res_function = resolution_Gauss;
% res_function=@(t,tau) resolution_Laplace_pleateau(t,tau,tau);

fitfunc="exp(-gamma*x/2)*(cos(omega*x)+gamma/omega/2*sin(omega*x))";

%% 0.5 Data assembly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Peak Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_period = 4;   % number of periods included in DO fit
t_cell=cell(n_runbin + 1,n_q - last_q_switch);
cf_cell=cell(n_runbin + 1,n_q - last_q_switch);
ft_cell=cell(n_runbin + 1,n_q - last_q_switch);
om_vals_cell=cell(n_runbin + 1,n_q - last_q_switch);
gamma_vals=zeros(n_runbin + 1,n_q - last_q_switch);
omega_1_vals=zeros(n_runbin + 1,n_q - last_q_switch);

fitob_cell=cell(n_runbin + 1,n_q - last_q_switch);
gamma_fit_vals=zeros(n_runbin + 1,n_q - last_q_switch);
omega_1_fit_vals=zeros(n_runbin + 1,n_q - last_q_switch);
gamma_confint_vals=zeros(2,n_runbin + 1,n_q - last_q_switch);
omega_1_confint_vals=zeros(2,n_runbin + 1,n_q - last_q_switch);
% gamma_confint_ste_vals=zeros(n_runbin + 1,n_q - last_q_switch);
% omega_1_confint_ste_vals=zeros(n_runbin + 1,n_q - last_q_switch);


ft_peak_vals=zeros(n_runbin + 1,n_q - last_q_switch);
om_peak_vals=zeros(n_runbin + 1,n_q - last_q_switch);
chi_vals_fit=zeros(n_runbin + 1,n_q - last_q_switch);

ft_spline_peak_vals=zeros(n_runbin + 1,n_q - last_q_switch);
ft_spline_om_vals=zeros(n_runbin + 1,n_q - last_q_switch);
ft_hwhm_vals=zeros(n_runbin + 1,n_q - last_q_switch);

for i_q=1:n_q-last_q_switch
    if isfile(fullfile_betterRes) && q_vals(i_q) > q_switch
        t=S_betterres.averaging_times;
        n_t=numel(t);
        q_indices=i_q:n_q:n_t*n_q;
        cf = cf_full_betterres(q_indices);
        cf_runs = real(cf_collect_betterres(:,q_indices));
    else
        t=S.averaging_times;
        n_t=numel(t);
        q_indices=i_q:n_q:n_t*n_q;
        cf = cf_full(q_indices);
        cf_runs = real(gmperpmperp_collect(:,q_indices))/sqrtN^2;
    end
    
    c = fit_DampedOscillator_RealSpace(t,cf/cf(1),n_period,weightexp,'omega_1');
    
    gamma_vals(1,i_q)=c(1);
    omega_1_vals(1,i_q)=c(2);
    chi_vals_fit(1,i_q)=cf(1);

    fitob = fit(t(:),cf(:)/cf(1),fitfunc,'StartPoint',[gamma_vals(1,i_q),omega_1_vals(1,i_q)]);
    fitob_cell{1,i_q} = fitob;
    gamma_fit_vals(1,i_q)=fitob.gamma;
    omega_1_fit_vals(1,i_q)=fitob.omega;
    ci=confint(fitob,.95);
    gamma_confint_vals(:,1,i_q)=ci(:,1);
    omega_1_confint_vals(:,1,i_q)=ci(:,2);


    tau=res_factor/c(1);
    t_cf=t;
    while tau > max(t_cf)
        dt=(t_cf(end)-t_cf(1))/(numel(t_cf)-1);
        [t_cf,cf] = continue_cf(t_cf,cf/cf(1),2,max(3*dt,omega_1_vals(i_q)/8/pi));
        cf = cf * chi_vals_fit(1,i_q);
    end
    t_cell{1,i_q}=t_cf;
    cf_cell{1,i_q}=cf;
    res_vals=res_function(t_cf,tau);

    [ft_vals,om_vals]=FT_correlation(t_cf, (cf .* res_vals), 0);
    ft_vals=real(ft_vals);
    ft_cell{1,i_q}=ft_vals;
    om_vals_cell{1,i_q}=om_vals;
    [ft_max,i_max]=max(ft_vals);
    om_peak_vals(1,i_q)=abs(om_vals(i_max));
    ft_peak_vals(1,i_q)=ft_max;
    
%         ft_spline=spline(om_vals(:),ft_vals(:));
%     ft_spline_halfheight=@(x) fnval(ft_spline,x) - ft_max/2;
    
    ft_spline_neg=spline(om_vals(:),-ft_vals(:));
    [splinemax,splinemaxarg]=fnmin(ft_spline_neg);
    ft_spline_peak_vals(1,i_q)=abs(splinemax);
    ft_spline_om_vals(1,i_q)=abs(splinemaxarg);
    ft_spline_halfheight=spline(om_vals(:),ft_vals(:)-abs(splinemax)/2);
    ft_spline_zeros=fnzeros(ft_spline_halfheight);
%         ft_hwhm_vals(1,i_q)=(ft_spline_zeros(1,4)-ft_spline_zeros(1,3))/2;
    ft_hwhm_vals(1,i_q)=max(abs(ft_spline_zeros(1,:)))-abs(splinemaxarg);


    for i_run = 1:n_runbin
%         t=S.averaging_times;
        run_indices=1+(i_run-1)*runsep:i_run*runsep;

%         cf = real(mean(gmperpmperp_collect(run_indices,q_indices)))/sqrtN^2;
        cf = mean(cf_runs(run_indices,:));
        c = fit_DampedOscillator_RealSpace(t,cf/cf(1),n_period,weightexp,'omega_1');
    
        gamma_vals(i_run + 1,i_q)=c(1);
        omega_1_vals(i_run + 1,i_q)=c(2);
        chi_vals_fit(i_run + 1,i_q)=cf(1);
    
        fitob = fit(t(:),cf(:)/cf(1),fitfunc,'StartPoint',[gamma_vals(i_run + 1,i_q),omega_1_vals(i_run + 1,i_q)]);

        fitob_cell{i_run + 1,i_q} = fitob;
        gamma_fit_vals(i_run + 1,i_q)=fitob.gamma;
        omega_1_fit_vals(i_run + 1,i_q)=fitob.omega;
        ci=confint(fitob,.95);
        gamma_confint_vals(:,i_run + 1,i_q)=ci(:,1);
        omega_1_confint_vals(:,i_run + 1,i_q)=ci(:,2);

        tau=res_factor/c(1);
        t_cf=t;
        while tau > max(t_cf)
            dt=(t_cf(end)-t_cf(1))/(numel(t_cf)-1);
            [t_cf,cf] = continue_cf(t_cf,cf/cf(1),2,max(3*dt,omega_1_vals(i_q)/8/pi));
            cf = cf * chi_vals_fit(1,i_q);
        end
        t_cell{i_run + 1,i_q}=t_cf;
        cf_cell{i_run + 1,i_q}=cf;
        res_vals=res_function(t_cf,tau);
    
        [ft_vals,om_vals]=FT_correlation(t_cf, (cf .* res_vals), 0);
        ft_vals=real(ft_vals);
        ft_cell{i_run + 1,i_q}=ft_vals;
        om_vals_cell{i_run + 1,i_q}=om_vals;
        [ft_max,i_max]=max(ft_vals);
        om_peak_vals(i_run + 1,i_q)=abs(om_vals(i_max));
        ft_peak_vals(i_run + 1,i_q)=ft_max;

        ft_spline_neg=spline(om_vals(:),-ft_vals(:));
        [splinemax,splinemaxarg]=fnmin(ft_spline_neg);
        ft_spline_peak_vals(i_run + 1,i_q)=abs(splinemax);
        ft_spline_om_vals(i_run + 1,i_q)=abs(splinemaxarg);
        ft_spline_halfheight=spline(om_vals(:),ft_vals(:)-abs(splinemax)/2);
        ft_spline_zeros=fnzeros(ft_spline_halfheight);
%         ft_hwhm_vals(i_run + 1,i_q)=(ft_spline_zeros(1,4)-ft_spline_zeros(1,3))/2;
        ft_hwhm_vals(i_run + 1,i_q)=max(abs(ft_spline_zeros(1,:)))-abs(splinemaxarg);


    end
end
gamma_confint_ste_vals=squeeze((gamma_confint_vals(2,:,:) - gamma_confint_vals(2,:,:))/3.92);
omega_1_confint_ste_vals=squeeze((omega_1_confint_vals(2,:,:) - omega_1_confint_vals(2,:,:))/3.92);

















%% 0.75 Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall settings, annotatiosn etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.8*columnwidth_cm]);






%% 1 Subplot top: omega fit vs FFT
ax_top=axes('InnerPosition',[.2 .5 .75 .4]);

% pos_plottop=get(gca,'Position');



n_period = 4;
om_max = .06;
linestyle='-';
linewidth=.5;
markersize=5;

colors={'r', 'b'};
markers={'s', 'd'};
dispnames={"$\omega_1$ fit",...
    "$\omega_{\textrm{peak}}$"};
% y_data={omega_1_vals(1,:),om_peak_vals(1,:)};
y_data={mean(omega_1_vals(2:end,:)),mean(om_peak_vals(2:end,:))};
y_errs={std(omega_1_vals(2:end,:))/sqrt(n_runbin),std(om_peak_vals(2:end,:))/sqrt(n_runbin)};
y_min = min(min(y_data{1},y_data{2}));
y_max = max(max(y_data{1},y_data{2}));
for i = 1:2
    y=y_data{i};
    y_err=y_errs{i};
%     h_plot{i} = errorbar(q_vals,y,y_err,...
%         'DisplayName',dispnames{i},...
%         'LineStyle', linestyle, ...
%         'LineWidth', linewidth, ...
%         'HandleVisibility', 'on', ...
%         'Color',colors{i});
    h_plot{i} = plot(q_vals,y,...
        'DisplayName',dispnames{i},...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'Marker', markers{i},...
        'MarkerSize', markersize,...
        'HandleVisibility', 'on', ...
        'Color',colors{i});
    hold on;
end
qq=linspace(.25,.4);
fit_pow_om=@(q) .15*q;

h_powerlaw_omega = plot(qq,fit_pow_om(qq),...
    'LineStyle', '--', ...
    'LineWidth', 1.2*linewidth, ...
    'HandleVisibility', 'off', ...
    'Color','k');
annotation_str = '$\propto q$';
h_qpower_annotation = text(qq(end/2), fit_pow_om(qq(end/2)),annotation_str,...
    'interpreter','latex',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top',...
    'fontsize',fontsize_labels,...
    'Color','black');

xlim([0 Inf])
% ylim([.015 .17]);
ylim([.8*y_min 1.5*y_max]);

hXLabel = xlabel('$q$','interpreter','latex');
hYLabel = ylabel('Frequency','interpreter','latex');
hLegend = legend('Location', 'South','interpreter','latex',...
    'NumColumns',1);
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'XGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YTick',0:.02:.14,...
    'YScale','log', 'XScale','log',...
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






annotation_str = {modeltag,sprintf('$N = (%d)^2$',sqrtN),sprintf('$T = %.2f$',T)};
pos = get(gca,'Position');
dim=[pos(1)+.02, ...
    pos(2)+pos(4)-.1, ...
    .125,.09];
h_annot=annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');



























%% 1 Subplot bottom: gamma vs peak height
ax_bot=axes('InnerPosition',[.2 .1 .75 .4]);
hold on;


n_period = 4;
om_max = .06;
linestyle='-';
linewidth=.5;
markersize=5;

colors={'r', 'b', '#0B0'};
markers={'s', 'd', 'o'};
dispnames={"$\gamma(q)$",...
    "$2\chi_{m\perp}(q)/S_{m\perp}^{\textrm{max}}(q)$",...
    "FWHM"};
y_data={gamma_vals(1,:),2./ft_peak_vals(1,:).*chi_init};
y_data={mean(gamma_vals(2:end,:)),...
    2./mean(ft_peak_vals(2:end,:)./chi_vals_fit(2:end,:)),...
    2*mean(ft_hwhm_vals(2:end,:))};
y_errs={std(gamma_vals(2:end,:))/sqrt(n_runbin),...
    std(2./ft_peak_vals(2:end,:).*chi_vals_fit(2:end,:))/sqrt(n_runbin),...
    2*std(ft_hwhm_vals(2:end,:))/sqrt(n_runbin)};
y_min = min([y_data{1},y_data{2},y_data{3}]);
y_max = max([y_data{1},y_data{2},y_data{3}]);
for i = 1:3
    y=y_data{i};
    y_err=y_errs{i};
    h_plot{i} = errorbar(q_vals,y,y_err,...
        'DisplayName',dispnames{i},...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'HandleVisibility', 'on', ...
        'Color',colors{i});
%     h_plot{i} = plot(q_vals,y,...
%         'DisplayName',dispnames{i},...
%         'LineStyle', linestyle, ...
%         'LineWidth', linewidth, ...
%         'Marker', markers{i},...
%         'MarkerSize', markersize,...
%         'HandleVisibility', 'on', ...
%         'Color',colors{i});
    fitobs{i} = fit(q_vals(:),y(:),'a*x^b');
end
xlim([0 Inf]);
% ylim([1e-3 1e-1]);
ylim([.8*y_min 1.5*y_max]);
% ylim([2e1 1e6]);
y=y_data{2};
fit_square = fit(q_vals(:),y(:),'a*x^2');
qq=linspace(.25,.4);
% fit_pow_gamma=@(q) .5*fitobs{1}.a*q.^(fitobs{1}.b);
fit_pow_gamma=@(q) .5*fit_square.a*q.^2;
% fit_pow_gamma=@(q) .5*.15*q.^(2);

h_powerlaw_gamma = plot(qq,fit_pow_gamma(qq),...
    'LineStyle', '--', ...
    'LineWidth', 1.2*linewidth, ...
    'HandleVisibility', 'off', ...
    'Color','k');
annotation_str = '$\propto q^2$';
h_qpower_annotation = text(qq(end/2), fit_pow_gamma(qq(end/2)),annotation_str,...
    'interpreter','latex',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top',...
    'fontsize',fontsize_labels,...
    'Color','black');

hXLabel = xlabel('$q$','interpreter','latex');
% hYLabel = ylabel('$S_{m\perp}(q,\omega) / \chi_{m\perp}(q)$','interpreter','latex');
% hYLabel = ylabel('$\gamma,S^{\textrm{max}}$','interpreter','latex');
hYLabel = ylabel('Damping','interpreter','latex');
hLegend = legend('Location', 'South','interpreter','latex',...
    'NumColumns',1);
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'XGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','log', 'XScale','log',...
    'XAxisLocation','bottom',...
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











    


%% 1 Saving
figname=sprintf('%s/%s_omegagamma_Compare_Fits_sqrtN_%d_T_%.3f',basedir,model,sqrtN,T);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end
















%% 2 Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall settings, annotatiosn etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.8*columnwidth_cm]);






%% 2 Subplot top: omega fit vs FFT
ax_top=axes('InnerPosition',[.2 .5 .75 .4]);

% pos_plottop=get(gca,'Position');



n_period = 4;
om_max = .06;
linestyle='-';
linewidth=.5;
markersize=5;

colors={'r', 'b'};
markers={'s', 'd'};
dispnames={"$\omega_1$ fit",...
    "$\omega_{\textrm{peak}}$"};
% y_data={omega_1_vals(1,:),om_peak_vals(1,:)};
y_data={mean(omega_1_vals(2:end,:)),mean(om_peak_vals(2:end,:))};
y_errs={std(omega_1_vals(2:end,:))/sqrt(n_runbin),std(om_peak_vals(2:end,:))/sqrt(n_runbin)};

y_min = min(min(y_data{1},y_data{2}));
y_max = max(max(y_data{1},y_data{2}));
for i = 1:2
    y=y_data{i};
    y_err = y_errs{i};
    h_plot{i} = errorbar(q_vals,y./q_vals,y_err./q_vals,...
        'DisplayName',dispnames{i},...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'HandleVisibility', 'on', ...
        'Color',colors{i});
%     h_plot{i} = plot(q_vals,y./q_vals,...
%         'DisplayName',dispnames{i},...
%         'LineStyle', linestyle, ...
%         'LineWidth', linewidth, ...
%         'Marker', markers{i},...
%         'MarkerSize', markersize,...
%         'HandleVisibility', 'on', ...
%         'Color',colors{i});
    hold on;
end

xlim([0 Inf])
% ylim([.015 .17]);

hXLabel = xlabel('$q$','interpreter','latex');
hYLabel = ylabel('$\omega / q$','interpreter','latex');
hLegend = legend('Location', 'South','interpreter','latex',...
    'NumColumns',1);
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'XGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','log',...
    'XAxisLocation','top',...
    'LineWidth', .5)
% ylim([1e1 Inf]);
%     'YTick',0:.02:.14,...


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






annotation_str = {modeltag,sprintf('$N = (%d)^2$',sqrtN),sprintf('$T = %.2f$',T)};
pos = get(gca,'Position');
dim=[pos(1)+.02, ...
    pos(2)+pos(4)-.1, ...
    .125,.09];
h_annot=annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');



























%% 2 Subplot bottom: gamma vs peak height
ax_bot=axes('InnerPosition',[.2 .1 .75 .4]);
hold on;


n_period = 4;
om_max = .06;
linestyle='-';
linewidth=.5;
markersize=5;

colors={'r', 'b', '#0B0'};
markers={'s', 'd', 'o'};
dispnames={"$\gamma(q)$",...
    "$2\chi_{m\perp}(q)/S_{m\perp}^{\textrm{max}}(q)$",...
    "FWHM"};
y_data={gamma_vals(1,:),2./ft_peak_vals(1,:).*chi_init};
y_data={mean(gamma_vals(2:end,:)),...
    2./mean(ft_peak_vals(2:end,:)./chi_vals_fit(2:end,:)),...
    2*mean(ft_hwhm_vals(2:end,:))};
y_errs={std(gamma_vals(2:end,:))/sqrt(n_runbin),...
    std(2./ft_peak_vals(2:end,:).*chi_vals_fit(2:end,:))/sqrt(n_runbin),...
    2*std(ft_hwhm_vals(2:end,:))/sqrt(n_runbin)};
y_min = min([y_data{1},y_data{2},y_data{3}]);
y_max = max([y_data{1},y_data{2},y_data{3}]);
for i = 1:3
    y=y_data{i};
    y_err = y_errs{i};
    h_plot{i} = errorbar(q_vals,y./q_vals.^2,y_err./q_vals.^2,...
        'DisplayName',dispnames{i},...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'HandleVisibility', 'on', ...
        'Color',colors{i});
%     h_plot{i} = plot(q_vals,y./q_vals.^2,...
%         'DisplayName',dispnames{i},...
%         'LineStyle', linestyle, ...
%         'LineWidth', linewidth, ...
%         'Marker', markers{i},...
%         'MarkerSize', markersize,...
%         'HandleVisibility', 'on', ...
%         'Color',colors{i});
    fitobs{i} = fit(q_vals(:),y(:),'a*x^b');
end
xlim([0 Inf]);
% ylim([1e-3 1e-1]);
% ylim([.8*y_min 1.5*y_max]);
% ylim([2e1 1e6]);


hXLabel = xlabel('$q$','interpreter','latex');
% hYLabel = ylabel('$S_{m\perp}(q,\omega) / \chi_{m\perp}(q)$','interpreter','latex');
% hYLabel = ylabel('$\gamma,S^{\textrm{max}}$','interpreter','latex');
hYLabel = ylabel('$\gamma / q^2$','interpreter','latex');
hLegend = legend('Location', 'Best','interpreter','latex',...
    'NumColumns',1);
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'XGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','log',...
    'XAxisLocation','bottom',...
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
figname=sprintf('%s/%s_omegagamma_Scaled_Fits_sqrtN_%d_T_%.3f',basedir,model,sqrtN,T);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end