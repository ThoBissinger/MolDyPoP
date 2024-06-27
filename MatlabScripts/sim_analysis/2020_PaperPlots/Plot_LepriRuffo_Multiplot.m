%% 0 Initialize
run initialization_script;
basedir=sprintf('%s/plots/LepriRuffo',fig_base);

saveswitch = 1;
    

%% 0 Preparing data
% model = "mxy"; modeldir = "mxy_3.00";
% mxy_dirs= "/data/scc/thobi/211201_LongerTime/mxy_3.00";
modelnames=["mxy" "fmxy"];
modelnames_text=["MXY" "DXY"];
data_rootdir="/data/scc/thobi";
data_dirs={"220201_ReducedSmapleStepDeltat", "210727_LepriRuffo_GoodResolution",  "211201_LongerTime", ...
    "210715_LinearTimeSampling"};
sampfilenames={"samp_Dynamics", "samp_LepriRuffo", "samp_Dynamics", "samp_Dynamics"};
data_dirs={"210727_LepriRuffo_GoodResolution", "220201_ReducedSmapleStepDeltat", "211201_LongerTime", ...
    "210715_LinearTimeSampling"};
sampfilenames={"samp_LepriRuffo", "samp_Dynamics", "samp_Dynamics", "samp_Dynamics"};

% fmxy_dir= "/data/scc/thobi/211201_LongerTime/fmxy";

sqrtN_vals = [16 32 64 128 256];
L_vals = [9.25 18.5 37 74 148];
N_N = numel(sqrtN_vals);

T_dirs = {"T_.03" "T_.17"};
T_vals = [.03 .17];
offset_scale = [2 2];
N_T = numel(T_vals);

i_T_short = 2; % index to have the zoomed in short-time version
labels = ["$(a)$", "$(b)$"; "$(c)$", "$(d)$"];

mxy_runmax=500;
fmxy_runmax=500;

y_min_vals = [.9, .5];
y_max_vals = [1.1, 1.5];

n_figspacing = 3;
x_rescaling_ratio = .9;
y_rescaling_ratio = .8;

% 
% marker_types=["v" "^" "o" "d"];
plots_xwidth = .38;
plots_ywidth = .32;
plots_xpos = .5 - plots_xwidth * [1, 0];
plots_ypos = flip(.55 - plots_ywidth * [1, 0]); 
% plots_ypos = .9 - plots_ywidth*(1:2); 

data_absM=zeros(2,N_N,N_T);
data_t=cell(2,N_N,N_T);
data_cf=cell(2,N_N,N_T);
data_eta_vals = zeros(2,N_T);
for i_N = 1:N_N
    sqrtN=sqrtN_vals(i_N);
    N = sqrtN^2;
    for i_T = 1:N_T
%         curfile=sprintf('%s/sqrtN_%d/%s/samp_Dynamics.mat',mxy_dir,sqrtN,T_dirs{i_T});
        i_dir=1;
        curfile=sprintf('%s/%s/mxy_3.00/sqrtN_%d/%s/%s.mat',data_rootdir,data_dirs{i_dir},sqrtN,T_dirs{i_T},sampfilenames{i_dir});
        while (~ isfile(curfile))
            i_dir = i_dir + 1;
            curfile=sprintf('%s/%s/mxy_3.00/sqrtN_%d/%s/%s.mat',data_rootdir,data_dirs{i_dir},sqrtN,T_dirs{i_T},sampfilenames{i_dir});
        end
        load(curfile,'averaging_times','absM_av','ACF_Spin');
        data_absM(1,i_N,i_T) = absM_av;
        data_t{1,i_N,i_T} = averaging_times;
        data_cf{1,i_N,i_T} = ACF_Spin;
        
        i_dir=1;
        curfile=sprintf('%s/%s/fmxy/sqrtN_%d/%s/%s.mat',data_rootdir,data_dirs{i_dir},sqrtN,T_dirs{i_T},sampfilenames{i_dir});
        while (~ isfile(curfile))
            i_dir = i_dir + 1;
            curfile=sprintf('%s/%s/fmxy/sqrtN_%d/%s/%s.mat',data_rootdir,data_dirs{i_dir},sqrtN,T_dirs{i_T},sampfilenames{i_dir});
        end
        load(curfile,'averaging_times','absM_av','ACF_Spin');
        data_absM(2,i_N,i_T) = absM_av;
        data_t{2,i_N,i_T} = averaging_times;
        data_cf{2,i_N,i_T} = ACF_Spin;
    end
end
for i_T = 1:N_T
    eta_fitob=fit_eta_Magnetization_FS(data_absM(1,:,i_T),L_vals);
    data_eta_vals(1,i_T) = eta_fitob.eta;

    eta_fitob=fit_eta_Magnetization_FS(data_absM(2,:,i_T),L_vals);
    data_eta_vals(2,i_T) = eta_fitob.eta;
end








%% 1 Subplot loop
fig = figure;
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.35*columnwidth_cm]);

c_map=colormap_sqrtN();
for i_T = 1:2
    T=T_vals(i_T);
    for i_model = 1:2
        i_fig = i_model + 2*i_T - 2;
%         i_figstart = i_model + 2*n_figspacing * (i_T - 1);
%         i_figend = i_figstart + 2*n_figspacing - 2;
%         subplot(2*n_figspacing+1,2,[i_figstart:2:i_figend],'replace');
        h_plot{i_model,i_T}=subplot(2,2,i_fig,'replace');

        for i_N = 1 : N_N
            sqrtN=sqrtN_vals(i_N);
            L = L_vals(i_N);
            
            dispname=sprintf('$N = (%d)^2$', sqrtN);
            t = data_t{i_model,i_N,i_T};
            cf = data_cf{i_model,i_N,i_T};
            eta = data_eta_vals(i_model,i_T);
            
            h_sim_data(i_N) = plot(t/L,L^(eta)*cf,...
                'LineStyle', '-', ...
                'LineWidth', .7, ...
                'DisplayName', dispname, ...
                'Color',c_map(i_N,:));
            hold on;
        %         'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k', ...
        
        end
        xlim([0,9.99]);
        ylim([y_min_vals(i_T),y_max_vals(i_T)]);
        
        
        hXLabel = xlabel('$t/L$','interpreter','latex');
        hYLabel = ylabel('$L^{\eta}C_{m}^{\textrm{inc}}(t)$','interpreter','latex');
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
            'LineWidth', .5, ...
            'XTick', 0:2:10, ...
            'XScale', 'linear','YScale','lin');
        
        if i_model == 1
            set(gca,'YAxisLocation','left');
        else
            set(gca,'YAxisLocation','right');
        end

        if i_T == 1
            set(gca,'XAxisLocation','top','XLabel',[]);
            set(gca,'YTick', .5:.05:3);
        else
            set(gca,'XAxisLocation','bottom');
            set(gca,'YTick', .4:.2:3);
        end

        pos=get(gca,'position');                            % retrieve the current values
        pos(1)=pos(1)+.5*(1-x_rescaling_ratio)*pos(3);      % recenters x-direction
        pos(3)=x_rescaling_ratio*pos(3);                    % scales width
        pos(2)=pos(2)+.5*(1-y_rescaling_ratio)*pos(4);      % recenters y-direction
        pos(4)=y_rescaling_ratio*pos(4);                    % scales height
        set(gca,'position',pos);                            % set to new values
    
        h_text=add_subfig_label(gca,labels(i_model,i_T),"se","lin","lin",fontsize_subfiglabels);
        

    end
end
    
for i_T = 1:2
    T=T_vals(i_T);
   
    for i_model = 1:2
        curpos = [plots_xpos(i_model), plots_ypos(i_T), plots_xwidth,plots_ywidth];
        set(h_plot{i_model,i_T},'Position',curpos);
        eta = data_eta_vals(i_model,i_T);

         str = {sprintf('$T = %.2f$',T),sprintf('$\\eta = %.2f$',eta)};
    %     pos=get(gca,'position');
        w_block = .15; w_block_spacing = .01;
        h_block = .09; h_block_spacing = .01;
        dim(1) = curpos(1) + curpos(3) - w_block - w_block_spacing;
        dim(2) = curpos(2) + curpos(4) - h_block - h_block_spacing;
        dim(3) = w_block;
        dim(4) = h_block;
    %         dim=[xpos_label(i_x) ypos_label(i_y) .055 .05];
        annotation('textbox',dim, 'String', str, 'FitBoxToText','off',...
            'interpreter','latex',...
            'LineWidth', .5, ...
            'VerticalAlignment','middle', 'HorizontalAlignment','left',...
            'Color','black','FontSize', fontsize_annotation,...
            'BackgroundColor','white');
        
    end
end




%% Legend subplot
% subplot(2*n_figspacing+1,2,[4*n_figspacing+1,4*n_figspacing+2],'replace');
% subplot(3,2,[5:6],'replace');
c_map=colormap_sqrtN();

pos=get(gca,'position');  % retrieve the current values
pos(3)=0.9*pos(3);        % try reducing height 10%
pos(4)=0.9*pos(4);        % try reducing width 10%
% set(gca,'position',pos);  % write the new values
% axis('off')
pow_legend = [0 0.03 1 .07];
hLegend = legend('Location', pow_legend,'interpreter','latex',...
    'NumColumns',1);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
























%% Overall labels and size setting


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall settings, annotatiosn etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_m_title=annotation('textbox',[.24 .93 .2 .05],...
    'string','MXY','Units','normalized',...
    'FitBoxToText','on',...
    'Interpreter','latex',...
    'EdgeColor','none',...
    'FontSize',fontsize_labels);

h_mperp_title=annotation('textbox',[.64 .93 .2 .05],...
    'string','DXY','Units','normalized',...
    'FitBoxToText','on',...
    'Interpreter','latex',...
    'EdgeColor','none',...
    'FontSize',fontsize_labels);

hLegend = legend('interpreter','latex',...
    'NumColumns',3);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
% hLegend.Position(1) = 0.05;
% hLegend.Position(2) = 0.08;


% set_fonts_default;























% Saving
if(saveswitch == 1)
    figname=sprintf('%s/ACF_Multiplot',basedir);
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end











%% 2 Short time plot
i_T = i_T_short;
% for i_T = 1:N_T
T = T_vals(i_T);
h_plot_st=cell(2,1);
fig=figure;
set(fig,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.6*columnwidth_cm]);
plot_xwidth = 0.7750;
plot_xpos = .13;
plot_ywidth = .42;
plot_ypos = [.08 .5];
pos_up = [plot_xpos,plot_ypos(2),plot_xwidth,plot_ywidth];
pos_down = [plot_xpos,plot_ypos(1),plot_xwidth,plot_ywidth];
for i_model = 1:2
    h_plot_st{i_model}=subplot(2,1,i_model,'replace');
    

    for i_N = 1 : N_N
        sqrtN=sqrtN_vals(i_N);
        L = L_vals(i_N);
        
        dispname=sprintf('$N = (%d)^2$', sqrtN);
    
        t = data_t{i_model,i_N,i_T};
        cf = data_cf{i_model,i_N,i_T};
        eta = data_eta_vals(i_model,i_T);
    
        h_sim_data(i_N) = plot(t/L,L^(eta)*cf,...
            'LineStyle', '-', ...
            'LineWidth', .5, ...
            'DisplayName', dispname, ...
            'Color',c_map(i_N,:));
        hold on;
    
    end
    offset_scale = [1.1 1.1];
    xx=linspace(.2,1.0);
    dispname='$\sim (t/L)^{-\eta}$';
    h_powerlaw=plot(xx,offset_scale(i_T)*xx.^-eta,'--','DisplayName',dispname,...
        'Color','black','LineWidth',2);
    xlim([4.4e-3,2]);
%     ylim([y_min_vals(i_T),y_max_vals(i_T)]);
    
    
    hXLabel = xlabel('$t/L$','interpreter','latex');
    hYLabel = ylabel('$L^{\eta}C_{m}^{\textrm{inc}}(t)$','interpreter','latex');
    h_axis = gca;
%     if i_model == 1 
%         hLegend = legend('Location', 'Southwest','interpreter','latex','FontSize', 20,...
%             'NumColumns',1);
%         set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
%     end
    % Font
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    % set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
    
    % Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'XGrid', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'LineWidth', .5,...
        'XScale', 'log','YScale','log');

%     str = {sprintf('$T = %.2f$',T),sprintf('$\\eta = %.2f$',eta)};
    pos=get(gca,'position');
    dim(1) = pos(1) + pos(3)/2;
%     dim(2) = pos(2) + pos(4);
    dim(1) = .5;
    dim(2) = .9;
    dim(3) = .02;
    dim(4) = .08;
    str=modelnames_text(i_model);
%         dim=[xpos_label(i_x) ypos_label(i_y) .055 .05];
    h_annot{i_model}=text('Position',[dim(1),dim(2)], 'String', str, ...
        'Units','normalized',...
        'FontWeight','bold',...
        'interpreter','latex',...
        'LineWidth', .5, ...
        'VerticalAlignment','middle', 'HorizontalAlignment','center',...
        'Color','black','FontSize', fontsize_titles,...
        'BackgroundColor','white');

    h_text=add_subfig_label(gca,labels(1,i_model),"ne","log","log",fontsize_subfiglabels);

    if i_model == 1
        hLegend = legend('Location', 'Southwest','interpreter','latex','FontSize', 12,...
            'NumColumns',1);
        set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

        set(gca,'XAxisLocation','top','XLabel',[]);
    else
        set(gca,'XAxisLocation','bottom');
    end
end
% plot_xwidth = 0.7750;
% plot_xpos = .13;
% plot_ywidth = .42;
% plot_ypos = [.08 .5];

% pos_up = get(h_plot_st{1},'Position');
% pos_down = get(h_plot_st{2},'Position');
% pos_down(2) = pos_up(2)-pos_down(4);
% set(h_plot_st{2},'Position',pos_down);
set(h_plot_st{1},'Position',pos_up);
set(h_plot_st{2},'Position',pos_down);










% Saving
if(saveswitch == 1)
%     figname=sprintf('%s/%s_LepriRuffo_ShortTime_T_%.3f',basedir,curmodel,T);
    figname=sprintf('%s/LepriRuffo_ShortTime_Multiplot_T_%.3f',basedir,T);
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end





% end