%% 0 Initialization
clear
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/chi_plots/Multiplot',fig_base);

%% 0 Preparing data
model = "mxy"; modeldir = "mxy_3.00";
% model = "fmxy"; modeldir = "fmxy";
% model = "xy"; modeldir = "xy_s";

sqrtN_vals = [16 32 64 128 256];
L_vals = [9.25 18.5 37 74 148];

linewidth = .8;



if model == "xy"
    sqrtN_Tplot=128;
    L_Tplot = 128;
else
    sqrtN_Tplot=256;
    L_Tplot = L_vals(end);
    % T_vals = [.03 .07 .11 .14 .17 .18 .19 .20 .22];
    % T_dirs = {"T_.03" "T_.07" "T_.11" "T_.14" "T_.17" "T_.18" "T_.19" "T_.20" "T_.22"};
    T_vals = [.14 .17 .19 .22 .25];
    T_dirs = {"T_.14" "T_.17" "T_.19" "T_.22" "T_.25"};
    T_vals_select = [.17 .21];
    T_dirs_select = {"T_.17" "T_.21"};
    eta_vals=[.23 .25];
end
N_T = numel(T_vals);

innerpos_ymins = flip([.1 .35 .7]);
innerpos_ywidth = .25;
% innerpos_ymax  = [.35 .6 .95];
innerpos_xmins = [.1 .375 .65 ];
innerpos_xwidth = .275;
% innerpos_xmax  = [.375 .65 .925];



datadir=sprintf("/data/scc/thobi/211201_LongerTime/%s",modeldir);
datadir_alternate=sprintf("/data/scc/thobi/210715_LinearTimeSampling/%s",modeldir);
% datadir_alternate=sprintf("/data/scc/thobi/211201_LongerTime/%s",modeldir);

% fdatadir="/data/scc/thobi/211201_LongerTime/fmxy";
% fdatadir_shorttime="/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00";

mxy_runmax=500;
fmxy_runmax=500;

mxy_res_factor=3/4*log(mxy_runmax);
fmxy_res_factor=3/4*log(fmxy_runmax);
res_function = resolution_Gauss;

marker_types=["v" "^" "o" "d"];
cornames = ["\chi_{m}", "\chi_{m\parallel}", "\chi_{m\perp}"];

data_T_chimq=cell(3,N_T);
data_T_qbin=cell(1,N_T);
data_FS_chimq=cell(3,numel(sqrtN_vals),numel(T_vals_select));
data_FS_qbin=cell(numel(sqrtN_vals),numel(T_vals_select));
data_FS_absM=zeros(numel(sqrtN_vals),numel(T_vals_select));
% eta_vals=zeros(size(T_vals_select));
for i_T = 1:numel(T_vals)
%     mxy_curfile=sprintf('%s/sqrtN_%d/%s/samp_Dynamics.mat',datadir,sqrtN_Tplot,T_dirs{i_T});
%     if (~ isfile(mxy_curfile))
        mxy_curfile=sprintf('%s/sqrtN_%d/%s/samp_Dynamics.mat',datadir_alternate,sqrtN_Tplot,T_dirs{i_T});
%     end
    load(mxy_curfile,'qbin','chimxq_av','chimyq_av',...
        'chimparq_av','chimperpq_av');
    N = sqrtN_Tplot^2;
    data_T_chimq{1,i_T} = (chimxq_av + chimyq_av)/N;
    data_T_chimq{2,i_T} = chimparq_av/N;
    data_T_chimq{3,i_T} = chimperpq_av/N;
    data_T_qbin{i_T} = qbin;
end

for i_N = 1:numel(sqrtN_vals)
    sqrtN=sqrtN_vals(i_N);
    N = sqrtN^2;
    for i_T = 1:numel(T_vals_select)
        mxy_curfile=sprintf('%s/sqrtN_%d/%s/samp_Dynamics.mat',datadir_alternate,sqrtN,T_dirs_select{i_T});
        load(mxy_curfile,'qbin','chimxq_av','chimyq_av',...
            'chimparq_av','chimperpq_av','absM_av');
        data_FS_chimq{1,i_N,i_T} = (chimxq_av + chimyq_av)/N;
        data_FS_chimq{2,i_N,i_T} = chimparq_av/N;
        data_FS_chimq{3,i_N,i_T} = chimperpq_av/N;
        data_FS_qbin{i_N,i_T} = qbin;
        data_FS_absM(i_N,i_T) = absM_av;
    end
end
% for i_T = 1:numel(T_vals_select)
%     eta_fitob=fit_eta_Magnetization_FS(data_FS_absM(:,i_T),L_vals);
%     eta_vals(i_T) = eta_fitob.eta;
% end






































%% 0.5 Overall labels and size setting


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall settings, annotatiosn etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labels=["\chi_{m}(q)", "\chi_{m\parallel}(q)", "\chi_{m\perp}(q)"];
for i_cor = 1:3
    xskip=.04;
    yskip=.01;
    currentpos=[innerpos_xmins(i_cor)+innerpos_xwidth/2 - xskip,...
        (max(innerpos_ymins)+innerpos_ywidth+1)/2-yskip,...
        2*xskip,...
        2*yskip];
    curlabel=sprintf('$%s$',labels(i_cor));
    h_m_title{i_cor}=annotation('textbox',currentpos,...
        'string',curlabel,'Units','normalized',...
        'Interpreter','latex',...
        'VerticalAlignment','middle',...
        'HorizontalAlignment','center',...
        'FitBoxToText','on',...
        'EdgeColor','none',...
        'FontSize',fontsize_titles,...
        'FontWeight','bold');
%     currentpos=[innerpos_xmins(i_label)+innerpos_xwidth/2,...
%         (max(innerpos_ymins)+innerpos_ywidth+1)/2];
%     h_m_title{i_label}=text('Position',currentpos,...
%         'string',labels(i_label),'Units','normalized',...
%         'Interpreter','latex',...
%         'VerticalAlignment','middle',...
%         'HorizontalAlignment','center',...
%         'FontSize',fontsize_labels);
end

% h_m_title{i_label}=annotation('textbox',[.2   .92 .2 .05],...
%         'string','$\chi_{m}(q)$','Units','normalized',...
%         'FitBoxToText','on',...
%         'Interpreter','latex',...
%         'EdgeColor','none',...
%         'FontSize',fontsize_labels);
% h_mperp_title=annotation('textbox',[.48 .92 .2 .05],...
%     'string','$\chi_{m\perp}(q)$','Units','normalized',...
%     'FitBoxToText','on',...
%     'Interpreter','latex',...
%     'EdgeColor','none',...
%     'FontSize',fontsize_labels);
% 
% h_mpar_title=annotation('textbox',[.76 .92 .2 .05],...
%     'string','$\chi_{m\parallel}(q)$','Units','normalized',...
%     'FitBoxToText','on',...
%     'Interpreter','latex',...
%     'EdgeColor','none',...
%     'FontSize',fontsize_labels);

% set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
set(gcf,'units','centimeters','OuterPosition',[0 0 pagewidth_cm .95*pageheight_cm]);
% set_fonts_default;




























%% 1 Subplot 1: chim for different T
% y_data = cell{3,N_T};
subfig_labels = ["$(a)$" "$(b)$" "$(c)$"];
for i_cor = 1:3
    subplot(3,3,i_cor,'replace');
    inpos_cur = [innerpos_xmins(i_cor) innerpos_ymins(1) innerpos_xwidth innerpos_ywidth];
    set(gca,'InnerPosition',inpos_cur)
    c_map=turbo(N_T+1);c_map=c_map(2:end,:);
    for i_T = 1 : N_T
        T=T_vals(i_T);
        qbin_cur=data_T_qbin{i_T};
        chimq_cur = data_T_chimq{i_cor,i_T};
        
        dispname=sprintf('$T = %.3f$', T);
        
        h_sim_data(i_T) = plot(qbin_cur,chimq_cur,...
            'LineStyle', '-', ...
            'LineWidth', linewidth, ...
            'DisplayName', dispname, ...
            'Color',c_map(i_T,:));
        hold on;
    %         'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k', ...
    
    end
    h_axis = gca;
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
        
    xlim([0,max(qbin_cur)]);
    ylim([1e0,2e3]);
    

    % Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'XGrid', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'LineWidth', .5,...
        'XScale', 'linear','YScale','log');
    if i_cor == 1
        hLegend = legend('Location', 'SouthWest','interpreter','latex',...
            'NumColumns',1);
        set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
        hYLabel = ylabel('$\chi(q)$','interpreter','latex');
        set(hYLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

        annotation_str = {sprintf('$N = (%d)^2$',sqrtN_Tplot)};
        pos=get(gca,'Position');
        dim=[pos(1)+pos(3)-.1, pos(2)+pos(4)-.035, .06, .03];
    
        annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
            'interpreter','latex',...
            'LineWidth', .5, ...
            'VerticalAlignment','middle', 'HorizontalAlignment','left',...
            'Color','black','FontSize', fontsize_annotation,...
            'BackgroundColor','white');
        set(gca, 'YAxisLocation','left');
    elseif i_cor == 2
        set(gca, 'YTickLabel',[]);
    elseif i_cor == 3
        set(gca, 'YAxisLocation','right');
        hYLabel = ylabel('$\chi(q)$','interpreter','latex');
        set(hYLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

    end
    
    hXLabel = xlabel('$q$','interpreter','latex');
%     hYLabel = ylabel('$\chi(q)$','interpreter','latex');
    set(hXLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    
    
    h_text=add_subfig_label(gca,subfig_labels(i_cor),"nw","lin","log",fontsize_subfiglabels);
    
    
    
%     dim=[.15 .703 .2 .05];
%     annotation('textbox',dim, 'String', 'a', 'FitBoxToText','on',...
%         'interpreter','latex',...
%         'LineWidth', .5, ...
%         'VerticalAlignment','middle', 'HorizontalAlignment','left',...
%         'Color','black','FontSize', fontsize_annotation,...
%         'BackgroundColor','white');

end






%% 4 Subplot 4: chim FS at transition temperature
subfig_labels = ["$(d)$", "$(e)$", "$(f)$";...
    "$(g)$", "$(h)$", "$(i)$"];
for i_cor = 1:3
    for i_T = 1:2
        subplot(3,3,3+3*(i_T - 1) + i_cor,'replace');
        inpos_cur = [innerpos_xmins(i_cor) innerpos_ymins(1+i_T) innerpos_xwidth innerpos_ywidth];
        set(gca,'InnerPosition',inpos_cur)
    
        T=T_vals_select(i_T);
        eta=eta_vals(i_T);
        L_exponent = -(2-eta);
        N_N = numel(sqrtN_vals);
        c_map=colormap_sqrtN();
        for i_N = 1 : N_N
            sqrtN = sqrtN_vals(i_N);
            L = L_vals(i_N);
            
            qbin_cur=data_FS_qbin{i_N,i_T};
            chimq_cur = data_FS_chimq{i_cor,i_N,i_T};
            
            
            dispname=sprintf('$N = (%d)^2$', sqrtN);
            
            h_sim_data(i_N) = plot(qbin_cur*L,L^(L_exponent)*chimq_cur,...
                'LineStyle', '-', ...
                'LineWidth', linewidth, ...
                'Marker','^', ...
                'DisplayName', dispname, ...
                'Color',c_map(i_N,:));
            hold on;
        %         'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k', ...
        
        end
        % xlim([0,max(qbin_cur)]);
        xlim([6,52]);
        ylim([5e-4,3e-1]);

        set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'XGrid', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
            'LineWidth', .5,...
            'XScale', 'log','YScale','log');
        h_axis = gca;
        set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
        
        if i_cor == 1
            if i_T == 1
                hLegend = legend('Location', 'SouthWest','interpreter','latex',...
                    'NumColumns',1);
                set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
            end
            hYLabel = ylabel('$\chi(q)$','interpreter','latex');
            set(hYLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    
%             annotation_str = {sprintf('$N = (%d)^2$',sqrtN_Tplot)};
%             dim=[.245 .785 .2 .05];
%             dim=[.245 .88 .2 .05];
%             
%             
            pos=get(gca,'Position');
            dim=[pos(1)+pos(3)-.08, pos(2)+pos(4)-.035, .06, .03];
            annotation_str = {sprintf('$T = %.2f$',T),sprintf('$\\eta = %.2f$',eta)};
%             dim=[.14 .495 .2 .05];
            % dim=[.245 .88 .2 .05];
            annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                'interpreter','latex',...
                'LineWidth', .5, ...
                'VerticalAlignment','top', 'HorizontalAlignment','left',...
                'Color','black','FontSize', fontsize_annotation,...
                'BackgroundColor','white');
            set(gca, 'YAxisLocation','left');
            hYLabel = ylabel(sprintf('$L^{-(2-\\eta)}%s$',labels(i_cor)),'interpreter','latex');
            set(hYLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
        elseif i_cor == 2
            set(gca, 'YTickLabel',[]);
        elseif i_cor == 3
            set(gca, 'YAxisLocation','right');
            hYLabel = ylabel(sprintf('$L^{-(2-\\eta)}%s$',labels(i_cor)),'interpreter','latex');
            set(hYLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    
        end
        if i_T == 1
            set(gca, 'XAxisLocation','top');
            xlabel '';
        else
            hXLabel = xlabel('$qL$','interpreter','latex');
            set(hXLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
        end
%         hLegend = legend('Location', 'SouthWest','interpreter','latex',...
%             'NumColumns',1);
        
        
%         hYLabel = ylabel('$L^{-(2-\eta)}\chi_{m}(q)$','interpreter','latex');
        
        % Font
        
        
%         set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
        
        h_text_T{i_T,i_cor}=add_subfig_label(gca,subfig_labels(i_T,i_cor),"nw","log","log",fontsize_subfiglabels);
        % Adjust axes properties
%         set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
%             'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
%             'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
%             'LineWidth', .5,...
%             'XScale', 'log','YScale','log');
        
        
        
%         dim=[.302 .578 .2 .05];
%         annotation('textbox',dim, 'String', 'd', 'FitBoxToText','on',...
%             'interpreter','latex',...
%             'LineWidth', .5, ...
%             'VerticalAlignment','middle', 'HorizontalAlignment','left',...
%             'Color','black','FontSize', fontsize_annotation,...
%             'BackgroundColor','white');
        
        
        

    end
end


















%% Saving
if(saveswitch == 1)
    figname=sprintf('%s/%s_chi_Multiplot',basedir,model);
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end






%% Same plot but with different scaling of first row





























%% 0.5 Overall labels and size setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall settings, annotatiosn etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% innerpos_ymins = flip([.1 .375 .65]);
innerpos_ymins = flip([.08 .355 .64]);
innerpos_ywidth = .275;
labels=["\chi_{m}(q)", "\chi_{m\parallel}(q)", "\chi_{m\perp}(q)"];
for i_cor = 1:3
    xskip=.04;
    yskip=.01;
    currentpos=[innerpos_xmins(i_cor)+innerpos_xwidth/2 - xskip,...
        (max(innerpos_ymins)+innerpos_ywidth+1)/2-yskip,...
        2*xskip,...
        6*yskip];
    curlabel=sprintf('$%s$',labels(i_cor));
    h_m_title{i_cor}=annotation('textbox',currentpos,...
        'string',curlabel,'Units','normalized',...
        'Interpreter','latex',...
        'VerticalAlignment','middle',...
        'HorizontalAlignment','center',...
        'FitBoxToText','on',...
        'EdgeColor','none',...
        'FontSize',fontsize_titles,...
        'FontWeight','bold');
%     currentpos=[innerpos_xmins(i_label)+innerpos_xwidth/2,...
%         (max(innerpos_ymins)+innerpos_ywidth+1)/2];
%     h_m_title{i_label}=text('Position',currentpos,...
%         'string',labels(i_label),'Units','normalized',...
%         'Interpreter','latex',...
%         'VerticalAlignment','middle',...
%         'HorizontalAlignment','center',...
%         'FontSize',fontsize_labels);
end

% h_m_title{i_label}=annotation('textbox',[.2   .92 .2 .05],...
%         'string','$\chi_{m}(q)$','Units','normalized',...
%         'FitBoxToText','on',...
%         'Interpreter','latex',...
%         'EdgeColor','none',...
%         'FontSize',fontsize_labels);
% h_mperp_title=annotation('textbox',[.48 .92 .2 .05],...
%     'string','$\chi_{m\perp}(q)$','Units','normalized',...
%     'FitBoxToText','on',...
%     'Interpreter','latex',...
%     'EdgeColor','none',...
%     'FontSize',fontsize_labels);
% 
% h_mpar_title=annotation('textbox',[.76 .92 .2 .05],...
%     'string','$\chi_{m\parallel}(q)$','Units','normalized',...
%     'FitBoxToText','on',...
%     'Interpreter','latex',...
%     'EdgeColor','none',...
%     'FontSize',fontsize_labels);

% set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
set(gcf,'units','centimeters','OuterPosition',[0 0 .9*pagewidth_cm .8*pageheight_cm]);
% set_fonts_default;




























%% 1 Subplot 1: chim for different T
% y_data = cell{3,N_T};
subfig_labels = ["$(a)$" "$(b)$" "$(c)$"];
for i_cor = 1:3
    subplot(3,3,i_cor,'replace');
    inpos_cur = [innerpos_xmins(i_cor) innerpos_ymins(1) innerpos_xwidth innerpos_ywidth];
    set(gca,'InnerPosition',inpos_cur)
    c_map=turbo(N_T+1);c_map=c_map(2:end,:);
    for i_T = 1 : N_T
        T=T_vals(i_T);
        qbin_cur=data_T_qbin{i_T};
        chimq_cur = data_T_chimq{i_cor,i_T};
        
        dispname=sprintf('$T = %.3f$', T);
        
        h_sim_data(i_T) = plot(L_Tplot*qbin_cur,chimq_cur,...
            'LineStyle', '-', ...
            'LineWidth', linewidth, ...
            'DisplayName', dispname, ...
            'Color',c_map(i_T,:));
        hold on;
    %         'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k', ...
    
    end
    h_axis = gca;
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
        
%     xlim([0,max(qbin_cur)]);
    xlim([6,52]);
    ylim([7e-1,2e3]);
    

    % Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'XGrid', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'LineWidth', .5,...
        'XScale', 'log','YScale','log');
    if i_cor == 1
        hLegend = legend('Location', 'SouthWest','interpreter','latex',...
            'NumColumns',1);
        set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
        hYLabel = ylabel('$\chi(q)$','interpreter','latex');
        set(hYLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

        annotation_str = {sprintf('$N = (%d)^2$',sqrtN_Tplot)};
        pos=get(gca,'Position');
        dim=[pos(1)+pos(3)-.105, pos(2)+pos(4)-.035, .06, .03];
    
        annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
            'interpreter','latex',...
            'LineWidth', .5, ...
            'VerticalAlignment','middle', 'HorizontalAlignment','left',...
            'Color','black','FontSize', fontsize_annotation,...
            'BackgroundColor','white');
        set(gca, 'YAxisLocation','left');
    elseif i_cor == 2
        set(gca, 'YTickLabel',[]);
    elseif i_cor == 3
        set(gca, 'YAxisLocation','right');
        hYLabel = ylabel('$\chi(q)$','interpreter','latex');
        set(hYLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

    end
    set(gca, 'XAxisLocation','top');
    
    hXLabel = xlabel('$qL$','interpreter','latex');
%     hYLabel = ylabel('$\chi(q)$','interpreter','latex');
    set(hXLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    
    
    h_text=add_subfig_label(gca,subfig_labels(i_cor),"nw","log","log",fontsize_subfiglabels);
    
    
    
%     dim=[.15 .703 .2 .05];
%     annotation('textbox',dim, 'String', 'a', 'FitBoxToText','on',...
%         'interpreter','latex',...
%         'LineWidth', .5, ...
%         'VerticalAlignment','middle', 'HorizontalAlignment','left',...
%         'Color','black','FontSize', fontsize_annotation,...
%         'BackgroundColor','white');

end






%% 4 Subplot 4: chim FS at transition temperature
subfig_labels = ["$(d)$", "$(e)$", "$(f)$";...
    "$(g)$", "$(h)$", "$(i)$"];
for i_cor = 1:3
    for i_T = 1:2
        subplot(3,3,3+3*(i_T - 1) + i_cor,'replace');
        inpos_cur = [innerpos_xmins(i_cor) innerpos_ymins(1+i_T) innerpos_xwidth innerpos_ywidth];
        set(gca,'InnerPosition',inpos_cur)
    
        T=T_vals_select(i_T);
        eta=eta_vals(i_T);
        L_exponent = -(2-eta);
        N_N = numel(sqrtN_vals);
        c_map=colormap_sqrtN();
        for i_N = 1 : N_N
            sqrtN = sqrtN_vals(i_N);
            L = L_vals(i_N);
            
            qbin_cur=data_FS_qbin{i_N,i_T};
            chimq_cur = data_FS_chimq{i_cor,i_N,i_T};
            
            
            dispname=sprintf('$N = (%d)^2$', sqrtN);
            
            h_sim_data(i_N) = plot(qbin_cur*L,L^(L_exponent)*chimq_cur,...
                'LineStyle', '-', ...
                'LineWidth', linewidth, ...
                'Marker','^', ...
                'DisplayName', dispname, ...
                'Color',c_map(i_N,:));
            hold on;
        %         'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k', ...
        
        end
        % xlim([0,max(qbin_cur)]);
        xlim([6,52]);
        ylim([5e-4,3e-1]);

        set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'XGrid', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
            'LineWidth', .5,...
            'XScale', 'log','YScale','log');
        h_axis = gca;
        set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
        
        if i_cor == 1
            if i_T == 1
                hLegend = legend('Location', 'SouthWest','interpreter','latex',...
                    'NumColumns',1);
                set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
            end
            hYLabel = ylabel('$\chi(q)$','interpreter','latex');
            set(hYLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    
%             annotation_str = {sprintf('$N = (%d)^2$',sqrtN_Tplot)};
%             dim=[.245 .785 .2 .05];
%             dim=[.245 .88 .2 .05];
%             
%             
            pos=get(gca,'Position');
            dim=[pos(1)+pos(3)-.09, pos(2)+pos(4)-.035, .06, .03];
            annotation_str = {sprintf('$T = %.2f$',T),sprintf('$\\eta = %.2f$',eta)};
%             dim=[.14 .495 .2 .05];
            % dim=[.245 .88 .2 .05];
            annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                'interpreter','latex',...
                'LineWidth', .5, ...
                'VerticalAlignment','top', 'HorizontalAlignment','left',...
                'Color','black','FontSize', fontsize_annotation,...
                'BackgroundColor','white');
            set(gca, 'YAxisLocation','left');
%             hYLabel = ylabel(sprintf('$L^{-(2-\\eta)}%s$',labels(i_cor)),'interpreter','latex');
            hYLabel = ylabel('$L^{-(2-\eta)}\chi(q)$','interpreter','latex');
            set(hYLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
        elseif i_cor == 2
            set(gca, 'YTickLabel',[]);
        elseif i_cor == 3
            set(gca, 'YAxisLocation','right');
%             hYLabel = ylabel(sprintf('$L^{-(2-\\eta)}%s$',labels(i_cor)),'interpreter','latex');
            hYLabel = ylabel('$L^{-(2-\eta)}\chi(q)$','interpreter','latex');
            set(hYLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    
        end
        if i_T == 1
            set(gca, 'XAxisLocation','top');
            set(gca, 'XTickLabel',[]);
            xlabel '';
        else
            hXLabel = xlabel('$qL$','interpreter','latex');
            set(hXLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
        end
%         hLegend = legend('Location', 'SouthWest','interpreter','latex',...
%             'NumColumns',1);
        
        
%         hYLabel = ylabel('$L^{-(2-\eta)}\chi_{m}(q)$','interpreter','latex');
        
        % Font
        
        
%         set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
        
        h_text_T{i_T,i_cor}=add_subfig_label(gca,subfig_labels(i_T,i_cor),"nw","log","log",fontsize_subfiglabels);
        % Adjust axes properties
%         set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
%             'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
%             'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
%             'LineWidth', .5,...
%             'XScale', 'log','YScale','log');
        
        
        
%         dim=[.302 .578 .2 .05];
%         annotation('textbox',dim, 'String', 'd', 'FitBoxToText','on',...
%             'interpreter','latex',...
%             'LineWidth', .5, ...
%             'VerticalAlignment','middle', 'HorizontalAlignment','left',...
%             'Color','black','FontSize', fontsize_annotation,...
%             'BackgroundColor','white');
        
        
        

    end
end


















%% Saving
if(saveswitch == 1)
    figname=sprintf('%s/%s_chi_Multiplot',basedir,model);
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end


