clear;
close all;
initialization_script;

% fontsize_annotation = 8;
% fontsize_axis = 7;
% fontsize_ax_labels = 8;

basedir=sprintf('%s/plots/Magnetization_Evolution/EvolutionOnSphere',fig_base);


dirs=["/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_16",
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_64",
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256"];
sqrtN_vals = [16, 64, 256];
T_dirs=["T_.14", "T_.17", "T_.20"];
T_vals = [.14, .17, .20];
sampfilename = "sampling_output_Dynamics";
N_T = numel(T_vals);
N_N = numel(sqrtN_vals);
plots_xwidth = .25;
plots_ywidth = plots_xwidth;
plots_xpos = .5-1.5*plots_xwidth+(0:2)*plots_xwidth;
plots_ypos = flip(.5-1.5*plots_ywidth+(0:2)*plots_ywidth);


figure;
set(gcf,'units','centimeters','Position',[0 0 .8*pagewidth_cm .8*pagewidth_cm],...
    'Resize','off');

h_plots=cell(N_N,N_T);
for i_N = 1:3
    for i_T = 1:N_T
        
        curfile = sprintf("%s/%s/run_1/output/%s",dirs(i_N),T_dirs(i_T),sampfilename);
        sqrtN = sqrtN_vals(i_N);
        T = T_vals(i_T);
        load(curfile);
        t_max = max(averaging_times);
        n_t = numel(averaging_times);
        t_max_exp = round(log10(t_max));
        t_max_noexp = t_max / 10^t_max_exp;
        subfig_index = (i_N - 1) *N_T + i_T;
        h_plots{i_N,i_T}=subplot(N_N,N_T,subfig_index,'replace');
%         curpos = [plots_xpos(i_T),plots_ypos(i_N),plots_xwidth,plots_ywidth];
%         set(h_plots{i_N,i_T},'InnerPosition',curpos)
        hold on;
        plot(absM_av * cos(0:.01:2*pi),absM_av * sin(0:.01:2*pi),...
            'Color','red',...
            'LineWidth',2)
        Mx = zeros(1,3*n_t);
        My = zeros(1,3*n_t);
        for i_collect=1:3
            curfile = sprintf("%s/%s/run_%d/output/%s",dirs(i_N),T_dirs(i_T),(i_collect)*125 + 1, sampfilename);
            load(curfile,'M');
            Mx((i_collect - 1)*n_t + 1 : i_collect*n_t) = M(1,:);
            My((i_collect - 1)*n_t + 1 : i_collect*n_t) = M(2,:);
        end
        plot(Mx,My,...
            'Color','blue',...
            'LineWidth',.25); 
%         plot(M(1,:),M(2,:),...
%             'Color','blue',...
%             'LineWidth',.3); 
        
        ax_max = 1.2;
        xlim([-ax_max,ax_max]);
        ylim([-ax_max,ax_max]);
%         pbaspect([1 1 1])
    
    
        h_axis = gca;
        hXLabel = xlabel('$m_x$','interpreter','latex');
        hYLabel = ylabel('$m_y$','interpreter','latex');
    
        set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
        set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    
        set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'XGrid', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
            'XTick', -1:.5:1, 'YTick', -1:.5:1, ...
            'LineWidth', .5)

        annotation_str = sprintf('$\\langle m\\rangle = %.2f$',absM_av);
        h_mval = text(.9, .1,annotation_str,'units','normalized',...
            'interpreter','latex',...
            'HorizontalAlignment','right',...
            'fontsize',fontsize_annotation,...
            'Color','red');
%         annotation_str = {sprintf('$\\langle m\\rangle = %.2f$',absM_av)};
%         dim=[.56 .13 .1 .2];
%         annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
%             'interpreter','latex', 'Units','normalized',...
%             'VerticalAlignment','middle', 'HorizontalAlignment','left',...
%             'Color','black','FontSize', fontsize_annotation,...
%             'BackgroundColor','white');
        if (i_N == 1)
            h_Tlabel = text(.5, 1.3,sprintf("$T = %.2f$",T)','units','normalized',...
                'interpreter','latex',...
                'HorizontalAlignment','center',...
                'fontsize',fontsize_titles);
            set(gca,'XAxisLocation','top');
        elseif (i_N == 2)
            set(gca,'XLabel',[],'XTickLabel',[]);
        else
            set(gca,'XAxisLocation','bottom');
        end
        if (i_T == 1)
            h_Nlabel = text(-.4, .5,sprintf("$N = (%d)^2$",sqrtN)','units','normalized',...
                'interpreter','latex',...
                'HorizontalAlignment','center',...
                'fontsize',fontsize_titles,...
                'Rotation',90);
%             set(h_Nlabel,'Rotation',90);
            set(gca,'YAxisLocation','left');
        elseif (i_T == 2)
            set(gca,'YLabel',[],'YTickLabel',[]);
        else
            set(gca,'YAxisLocation','right');
        end
    end 
end

% Reassign size to avoid erasure during while creation
for i_N = 1:N_N
    for i_T = 1:N_T
        curpos = [plots_xpos(i_T),plots_ypos(i_N),plots_xwidth,plots_ywidth];
        set(h_plots{i_N,i_T},'InnerPosition',curpos)
    end
end
% gcf

% set(gcf,'units','centimeters','OuterPosition',[0 0 pagewidth_cm pageheight_cm]);

%     figname=sprintf('%s/mxy_sqrtN_%d_T_%.2f',basedir,sqrtN,T);
figname=sprintf('%s/mxy_M_Evolution',basedir);
fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

