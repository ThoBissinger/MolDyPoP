%% 0 Initialization
clear; close all;
run ../../initialization_script.m

S_shorttime = load('/data/scc/thobi/220125_SmallerTimeStep/mxy_3.00/sqrtN_128/T_.17/samp_LepriRuffo.mat');
S_longtime= load('/data/scc/thobi/210727_LepriRuffo_GoodResolution/mxy_3.00/sqrtN_128/T_.17/samp_LepriRuffo.mat');
S_verylongtime= load('/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_128/T_.17/samp_Dynamics.mat');

% t_s = S_shorttime.averaging_times(2:end);
% t_l = S_longtime.averaging_times(2:end);
% t_vl = S_verylongtime.averaging_times(2:end);
% 
% y_s = S_shorttime.ACF_MSD(2:end);
% y_l = S_longtime.ACF_MSD(2:end);
% y_vl = S_verylongtime.ACF_MSD(2:end);

t_s = S_shorttime.averaging_times+.001;
t_l = S_longtime.averaging_times+.01;
t_vl = S_verylongtime.averaging_times+.01;

y_s = S_shorttime.ACF_MSD;
y_l = S_longtime.ACF_MSD;
y_vl = S_verylongtime.ACF_MSD;

[t_l,y_l] = combine_cf(t_l,y_l,t_vl,y_vl);

%% 1 Combined plot
figure
% loglog(t,y,'-');
loglog(t_s,y_s,'-',...
    'Displayname','$dt = .001$',...
    'LineWidth',2.5);
hold on;
loglog(t_l,y_l,'-',...
    'Displayname','$dt = .01$',...
    'LineWidth',2.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust axes, legend properties etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend=legend('Location','southeast','Interpreter','latex',...
        'NumColumns',1);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', .5,'Xscale', 'log','Yscale', 'log')

hXLabel = xlabel('$t$','interpreter','latex');
hYLabel = ylabel('$\langle \Delta r^2(t)\rangle$','interpreter','latex');
h_axis = gca;

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

tt_short = logspace(-2, -.5,100);
y_fit_short = .3 * y_s(10)/t_s(10)^2 * tt_short.^2;
plot(tt_short,y_fit_short,'--',...
    'Color','black',...
    'LineWidth',2.5,...
    'HandleVisibility','off');
annotation_str = '$\propto t^{2}$';
h_tshort_annotation = text(tt_short(end/2), y_fit_short(end/2),annotation_str,...
    'interpreter','latex',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top',...
    'fontsize',fontsize_labels,...
    'Color','black');

tt_long = logspace(2, 4,100);
y_fit_long = .5 * y_l(end)/t_l(end) * tt_long;
plot(tt_long,y_fit_long,'--',...
    'Color','black',...
    'LineWidth',2.5,...
    'HandleVisibility','off');
annotation_str = '$\propto t$';
h_tlong_annotation = text(tt_long(end/2), y_fit_long(end/2),annotation_str,...
    'interpreter','latex',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top',...
    'fontsize',fontsize_labels,...
    'Color','black');

set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

% figname=sprintf('%s/mxy_MSD_ShortTime',basedir);
figname='mxy_MSD_ShortTime';

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end