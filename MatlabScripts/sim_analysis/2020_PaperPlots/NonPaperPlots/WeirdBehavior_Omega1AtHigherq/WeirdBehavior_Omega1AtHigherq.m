clear;
close all;
%% 0 Initialization
run ../../initialization_script.m
load('/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_128/T_.14/samp_Dynamics.mat',...
    'qbin','averaging_times','gmperpmperp');
t_highq=averaging_times;
g_highq=gmperpmperp;
q_vals_highq=qbin;
dq = qbin(1);
L = 2 * pi / dq;
qmax_highq = max(qbin);
n_q_highq=numel(qbin);

load('/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_128/T_.14/samp_Dynamics.mat',...
    'qbin','averaging_times','gmperpmperp');
t_lowq=averaging_times;
g_lowq=gmperpmperp;
q_vals_lowq=qbin;
qmax_lowq = max(qbin);
n_q_lowq=numel(qbin);






%% 1 Plot
n_period=4;
weightexp=1;
i_q = 3;
figure
cf=real(gmperpmperp(i_q:n_q:end));
c = fit_DampedOscillator_RealSpace(t,cf,n_period,weightexp,'omega_1');
gamma = c(1);
omega_1 = c(2);
plot(t,cf,'.',...
    'DisplayName','sim data');
hold on;
plot(t,fitfunc_DO(t,cf(1),c),'-',...
    'DisplayName','sim data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust axes, legend properties etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend=legend('Location','southeast','Interpreter','latex',...
    'NumColumns',1);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')

hXLabel = xlabel('$t$','interpreter','latex');
hYLabel = ylabel('$C_{m\perp}(q,t)$','interpreter','latex');
h_axis = gca;

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
