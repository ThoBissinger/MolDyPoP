% Makes a Multiplot of the structural properties of the system, including a
% plot of the direct correlation g(r), the MSD and two system snapshots at
% very low and moderately low temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
initialization_script;

saveswitch=1;
basedir=sprintf('%s/plots/structure',fig_base);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data collection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collection of g(r) for first plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = "mxy";
load /data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_64/T_.11/samp_Dynamics_gr.mat
sqrtN=256;
gr_basedir=sprintf('/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_%d',sqrtN);
gr_sampfilename='samp_Dynamics_gr';
MSD_basedir=sprintf('/data/scc/thobi/210727_LepriRuffo_GoodResolution/mxy_3.00/sqrtN_%d',sqrtN);
MSD_sampfilename='samp_LepriRuffo';
T_str=[ "T_.01"  "T_.03"  "T_.09"  "T_.15"  "T_.19"  "T_.31"];
T_vals=[.01 .03 .09 .15 .19 .31];

dt = .01;

snapshot_T_str=["T_.01" "T_.09"];
snapshot_T_vals=[.01 .09];
snapshot_basedir=sprintf('/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_%d',sqrtN);
snapshot_runnr=[3,1];
snapshot_sampfilename='snapshot_Dynamics_final';

% Directories for calculation of diffusion coefficient
% D_T_str = ["T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"];
% D_T_vals = [.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52];
% D_T_str = ["T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.16" "T_.17" "T_.18" "T_.19" "T_.20" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"];
% D_T_vals = [.01 .03 .05 .07 .09 .11 .13 .14 .15 .16 .17 .18 .19 .20 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52];
D_T_str = ["T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.15" "T_.17" "T_.19" "T_.21" "T_.23" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"];
D_T_vals = [.01 .03 .05 .07 .09 .11 .13 .15 .17 .19 .21 .23 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52];

D_basedir=sprintf('/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_%d',sqrtN);
D_sampfilename='samp_Dynamics';


N_T = numel(T_str);



x_border = .02;     % Border width in x direction (for positioning of objects)
y_border = .01;     % Border width in y direction (for positioning of objects)
sublabel_box = [.073 .025];  % Size of sublabel boxes (a, b, ...)

%% Data preparation
gr_cell=cell(1,N_T);
rbin_cell=cell(1,N_T);
MSD_cell=cell(1,N_T);
MSD_t_cell=cell(1,N_T);
for i_T = 1:N_T
    gr_curfile=sprintf('%s/%s/%s',gr_basedir,T_str{i_T},gr_sampfilename);
    load(gr_curfile);
    rbin_cell{i_T}=rbin;
    gr_cell{i_T}=gr;

    MSD_curfile=sprintf('%s/%s/%s',MSD_basedir,T_str{i_T},MSD_sampfilename);
    load(MSD_curfile,'averaging_times','ACF_MSD');
    MSD_t_cell{i_T}=averaging_times + dt;
    MSD_cell{i_T}=ACF_MSD;
end
D_vals = zeros(1,numel(D_T_vals));
for i_T = 1:numel(D_T_vals)
    MSD_curfile=sprintf('%s/%s/%s',D_basedir,D_T_str{i_T},D_sampfilename);
    load(MSD_curfile,'averaging_times','ACF_MSD');
    D_vals(i_T)=ACF_MSD(end)/4/(averaging_times(end) + dt); % factor 4 because it should be 2d
end
r_snap=cell(1,numel(snapshot_T_vals));
for i_T = 1:numel(snapshot_T_vals)
    curfile=sprintf('%s/%s/run_%d/output/%s.out',...
        snapshot_basedir,snapshot_T_str{i_T},snapshot_runnr(i_T),snapshot_sampfilename);
    [r_snap{i_T},~,~,~] = mxy_snapshot_extract(curfile,'r','mxy');
end













%% General settings of image
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 2.6*columnwidth_cm]); 

%% 1 Subplot 1: g(r) with low-T inset
c_map=turbo(N_T+1); c_map=c_map(2:end,:);
subplot(5,2,[1:4],'replace');
hold on;
dispname = sprintf('$T = %.2f$',T_vals(1));
plot(nan,nan,'-',...
    'LineStyle', '-', 'LineWidth', 1.2, ...
    'DisplayName', dispname, ...
    'Color', c_map(1,:));
for i_T = 2:N_T
    T = T_vals(i_T);
    dispname = sprintf('$T = %.2f$',T);
    plot(rbin_cell{i_T},gr_cell{i_T},'-',...
        'LineStyle', '-', 'LineWidth', 1.2, ...
        'DisplayName', dispname, ...
        'Color', c_map(i_T,:));
    
end
xlim([.35 1.8]);
ylim([0 2.1]);
hLegend=legend('Location','southeast','Interpreter','latex',...
    'NumColumns',2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust axes properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')

hXLabel = xlabel('$r$','interpreter','latex');
hYLabel = ylabel('$g(r)$','interpreter','latex');
h_axis = gca;

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Label a
% pos=get(gca,'InnerPosition');
% dim=[pos(1)+x_border pos(2)+pos(4)-y_border-sublabel_box(2) sublabel_box];
% annotation('textbox',dim, 'String', 'a', 'FitBoxToText','off',...
%     'interpreter','latex',...
%     'LineWidth', .5, ...
%     'VerticalAlignment','middle', 'HorizontalAlignment','left',...
%     'Color','black','FontSize', fontsize_annotation,...
%     'BackgroundColor','white');
% pos = get(gca,'Position');
% dim=[pos(1)+pos(3)-x_border-labelbox_size(1), ...
%     pos(2)+y_border, ...
%     labelbox_size(1),...
%     labelbox_size(2)];
% textpos=[.94 .06];
% textpos=[.06 .94];
% NW=[min(xlim) + .02*diff(xlim), max(ylim)-.02*diff(ylim)];
% text('Position',NW, 'String', '(a)', ...
%     'interpreter','latex',...
%     'HorizontalAlignment','left',...
%     'VerticalAlignment','top',...
%     'Color','black','FontSize', fontsize_subfiglabels);
h_text=add_subfig_label(gca,"$(a)$","nw","lin","lin",fontsize_subfiglabels);

% Annotation box for system size and density
pos=get(gca,'InnerPosition');
annotation_str = {sprintf('$N = (%d)^2$',sqrtN),'$\rho=3.00$'};
dim=[pos(1)+2*x_border pos(2)+y_border .15 2*sublabel_box(2)];
annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=get(gca,'Position');
inset_w = pos(3)*.5;                % Width of inset
inset_h = pos(4)*.32;               % Height of inset
pos = [pos(1) + pos(3) - inset_w - x_border,...     % Inset left edge
    pos(2) + pos(4) - inset_h - y_border,...        % Inset bottom edge
    inset_w,...                                     % Inset width
    inset_h];                                       % Inset height
ax_inset = axes('Position',pos);

set(ax_inset, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', .5, 'Xscale', 'lin')
maxgr=0;
for i_T = 1 : N_T
    T = T_vals(i_T);
    dispname = sprintf('$T = %.2f$',T);
    plot(rbin_cell{i_T},gr_cell{i_T},'-',...
        'LineStyle', '-', 'LineWidth', 1.2, ...
        'DisplayName', dispname, ...
        'Color', c_map(i_T,:));
    hold on;
    maxgr = max(maxgr,max(gr_cell{i_T}));
end
ylim([0 1.2*maxgr]);
xlim([0 4]);

NW = [min(xlim) max(ylim)]+[1.6*diff(xlim) -diff(ylim)]*0.02;
SE = [max(xlim) min(ylim)]+[-diff(xlim) .5*diff(ylim)]*0.02;
text(SE(1),SE(2),'$r$',...
    'VerticalAlignment','bottom','HorizontalAlignment','right',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_ax_labels);
text(NW(1),NW(2),'$g(r)$',...
    'VerticalAlignment','top','HorizontalAlignment','left',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_ax_labels);
set(gca, 'FontName', 'cmr12')
%     set([hXLabel, hYLabel], 'FontName', 'cmr12')
set([gca], 'FontSize', fontsize_axis)
    
























%% 2 Subplot 2: MSD with diffusion constant inset
subplot(5,2,[5:8],'replace');
c_map=turbo(N_T+1); c_map=c_map(2:end,:);
for i_T = 1:N_T
    T = T_vals(i_T);
    dispname = sprintf('$T = %.3g$',T);

    t = MSD_t_cell{i_T};
    y = MSD_cell{i_T};
    h_plot{i_T}=loglog(t,y,'-',...
        'LineWidth',1.2,...
        'Color',c_map(i_T,:),...
        'DisplayName',dispname);
    hold on;
end
xlim([2e-2 1e2]);
h_axis = gca;
%     h_title = title(sprintf("$T = %.3f$",T),'interpreter','latex');
hXLabel = xlabel('$t$','interpreter','latex');
hYLabel = ylabel('$\langle \Delta r^2(t)\rangle$','interpreter','latex');

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);


tt=linspace(1e-2,1e-1);
a = 1e-2;
y_fit = a*tt.^2;
h_plotshort = plot(tt,y_fit,'--',...
    'Color','black','LineWidth',1.2,...
    'HandleVisibility','off');
annotation_str = '$\propto t^{2}$';
h_tshort_annotation = text(tt(end/2), y_fit(end/2),annotation_str,...
    'interpreter','latex',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top',...
    'fontsize',fontsize_labels,...
    'Color','black');

tt=linspace(3e0,3e1);
D = 2e-2;
y_fit = D*tt;
h_plotlong = plot(tt,y_fit,'--',...
    'Color','black','LineWidth',1.2,...
    'HandleVisibility','off');
annotation_str = '$\propto t$';
h_tlong_annotation = text(tt(end/2), y_fit(end/2),annotation_str,...
    'interpreter','latex',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top',...
    'fontsize',fontsize_labels,...
    'Color','black');
    


% pos=get(gca,'Position');
% dim=[pos(1)+x_border pos(2)+pos(4)-y_border-sublabel_box(2) sublabel_box];
% annotation('textbox',dim, 'String', 'b', 'FitBoxToText','off',...
%     'interpreter','latex',...
%     'LineWidth', .5, ...
%     'VerticalAlignment','middle', 'HorizontalAlignment','left',...
%     'Color','black','FontSize', fontsize_annotation,...
%     'BackgroundColor','white');
% xlog = log(xlim);
% ylog = log(ylim);
% xratio = max(xlim)/min(xlim);
% yratio = max(ylim)/min(ylim);
% NW=[exp(min(xlog)+diff(xlog)*.02), exp(max(ylog)-diff(ylog)*.02)];
% text('Position',NW, 'String', '(b)', ...
%     'interpreter','latex',...
%     'HorizontalAlignment','left',...
%     'VerticalAlignment','top',...
%     'Color','black','FontSize', fontsize_subfiglabels);
h_text=add_subfig_label(gca,"$(b)$","nw","log","log",fontsize_subfiglabels);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=get(gca,'Position');
inset_w = pos(3)*.5;                % Width of inset
inset_h = pos(4)*.38;               % Height of inset
pos = [pos(1) + pos(3) - inset_w - x_border,...     % Inset left edge
    pos(2) + 2*y_border,...                         % Inset bottom edge
    inset_w,...                                     % Inset width
    inset_h];                                       % Inset height
ax_inset = axes('Position',pos);

set(ax_inset, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', .5, 'Xscale', 'lin')
maxgr=0;
plot(D_T_vals,D_vals,...
    'LineStyle', '-', 'LineWidth', 1, ...
    'Marker','s',...
    'MarkerSize',3,...
    'Color', 'black');
% ylim([0 1.2*maxgr]);
xlim([0 max(D_T_vals)]);

NW = [min(xlim) max(ylim)]+[1.6*diff(xlim) -diff(ylim)]*0.02;
SE = [max(xlim) min(ylim)]+[-diff(xlim) .5*diff(ylim)]*0.02;
text(SE(1),SE(2),'$T$',...
    'VerticalAlignment','bottom','HorizontalAlignment','right',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_ax_labels);
text(NW(1),NW(2),'$D$',...
    'VerticalAlignment','top','HorizontalAlignment','left',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_ax_labels);
set(gca, 'FontName', 'cmr12')
%     set([hXLabel, hYLabel], 'FontName', 'cmr12')
set([gca], 'FontSize', fontsize_axis)
    




%% 3 Subplot 3 & 4: Very low and slightly low temperature snapshots
abc_label=["(c)" "(d)"];
for i_T = 1:2
    T = snapshot_T_vals(i_T);
    subplot(5,2,8+i_T,'replace');
%     curfile=sprintf('%s/%s/run_%d/output/%s.out',...
%         snapshot_basedir,snapshot_T_str{i_T},snapshot_runnr(i_T),snapshot_sampfilename);
%     single_snap_fig(curfile,'mxy','pos');
    scatter(r_snap{i_T}(1,:),r_snap{i_T}(2,:),3,'filled');
%     plot(r_snap{i_T}(1,:),r_snap{i_T}(2,:),...
%         'LineStyle','none',...
%         'Marker');
%     pos_old=get(gca,'Position')
%     pbaspect([1 1 1]);
    curaspect=pbaspect;
    curaspect=curaspect/curaspect(1);
    xlim([0 12]);
    ylim([0 12*curaspect(2)]);
%     pbaspect([1 1 1]);

    % Label a
    pos=get(gca,'InnerPosition');
%     pos=get(gca,'InnerPosition');
    dim=[pos(1)+x_border, pos(2)+pos(4)-y_border-sublabel_box(2), sublabel_box];
    annotation('textbox',dim, 'String', abc_label(i_T), 'FitBoxToText','on',...
        'interpreter','latex',...
        'LineWidth', .5, ...
        'VerticalAlignment','middle', 'HorizontalAlignment','left',...
        'Color','black','FontSize', fontsize_annotation,...
        'BackgroundColor','white');
    
    % Annotation box for system size and density
    annotation_str = {sprintf('$T = %.2f$',T)};
    dim=[pos(1)+x_border pos(2)+y_border 4*sublabel_box(1) sublabel_box(2)];
    annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
        'interpreter','latex',...
        'LineWidth', .5, ...
        'VerticalAlignment','middle', 'HorizontalAlignment','left',...
        'Color','black','FontSize', fontsize_annotation,...
        'BackgroundColor','white');
    hXLabel = xlabel('$x$','interpreter','latex');
    hYLabel = ylabel('$y$','interpreter','latex');
    h_axis = gca;
    
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin',...
        'XTick',[0:2:100],'YTick',[0:2:100])

end










%% saving the data
if(saveswitch == 1)
    figname=sprintf('%s/mxy_strutcutre_sqrtN_%d',basedir,sqrtN);
    fprintf('Creating figure %s\n',figname)
    
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end
