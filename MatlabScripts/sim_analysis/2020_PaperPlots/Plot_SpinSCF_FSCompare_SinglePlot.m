%% Initialization
clear all
close all
initialization_script
basedir=sprintf('%s/plots/Spin_SCF',fig_base);
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;

% 
% mxydata=load('mxy/rho_3.00_eq.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat');
% xydata=load('xy/lf0_eq.mat'); xyfit=load('xy/lf0_CritExpFit.mat');

data_rootdir="/data/scc/thobi";
label='$(b)$';
labelboxdim=[.05 .065];
boxbound=[.015 .02];
for i_model = 1

    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        sqrtN_vals=[16, 32, 64, 128, 256];
        L_vals = 9.25/16*sqrtN_vals;
        symbols=['o','+','s','x','d'];
%         data_dirs={"211201_LongerTime/mxy_3.00", "210715_LinearTimeSampling/mxy_3.00"};
%         sampfilenames={"samp_Dynamics", "samp_Dynamics"};
        data_dirs={"210715_LinearTimeSampling/mxy_3.00"};
        sampfilenames={"samp_Dynamics"};
        T_vals=[.14 .17 .19];
        T_dirs=["T_.14" "T_.17" "T_.18"];
        eta_vals=[.15 .236 .32];
        offset_vals=[0 0 0];
        ylim_vals=[.55 2.5];

        xpos_T = [.8, .7, .8];
        ypos_T = [.28 .8 .66];%[.13+i_T*.24;];
        yshifts = [.11, -.30, -.11];
       
    else
        curmodel="xy";
        curtitle="SXY model";
        
        r_min = 3;
        r_max = 55;
    end
    N_T = numel(T_vals);
    N_N = numel(sqrtN_vals);

    %% 1 Gather data
    data_rbin=cell(N_N,N_T);
    data_cf=cell(N_N,N_T);
    data_absM=zeros(N_N,N_T);
    for i_N = 1:N_N
        sqrtN=sqrtN_vals(i_N);
        for i_T = 1:N_T
            T_dir=T_dirs(i_T);
            i_dir=1;
            curfile=sprintf('%s/%s/sqrtN_%d/%s/%s.mat',data_rootdir,data_dirs{i_dir},sqrtN,T_dir,sampfilenames{i_dir});
            while (~ isfile(curfile))
                i_dir = i_dir + 1;
                curfile=sprintf('%s/%s/sqrtN_%d/%s/%s.mat',data_rootdir,data_dirs{i_dir},sqrtN,T_dir,sampfilenames{i_dir});
            end
            load(curfile,'rbin','SCF_Spin_av','absM_av');
            data_rbin{i_N,i_T} = rbin;
            data_cf{i_N,i_T} = SCF_Spin_av;
            data_absM(i_N,i_T) = absM_av;
        end
    end
%     eta_vals=zeros(1,N_T);
%     for i_T = 1:N_T
%         eta_fitob=fit_eta_Magnetization_FS(data_absM(:,i_T),L_vals);
%         eta_vals(i_T)=eta_fitob.eta;
%     end

    %% Plotting
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    hold on
    c_map = turbo(10); c_map = c_map([2,5,9],:);
    for i_T = 1:N_T
        T = T_vals(i_T);
%         eta=eta_vals_FS(i_T);
        eta=eta_vals(i_T);
        for i_N = 1:N_N
            sqrtN = sqrtN_vals(i_N);
            r_vals = data_rbin{i_N,i_T};
            SCF = data_cf{i_N,i_T};
            nonzero_indices = find(SCF ~= 0);
            r_min = r_vals(2);
            r_max = r_vals(end);

            rr = linspace(r_vals(1),r_vals(end));

            dispname=sprintf('$N = (%d)^2$', sqrtN);
            offset = offset_vals(i_T);
            hSCF_line(i_N,i_T) = line(r_vals(nonzero_indices)/L_vals(i_N),SCF(nonzero_indices)*L_vals(i_N)^eta + offset);
            hold on;
            set(hSCF_line(i_N,i_T), ...
                'LineStyle', 'none', 'LineWidth', 1, ...
                'Marker', symbols(i_N), 'MarkerSize', 4, ...'Marker', '.', 'MarkerSize', 6, ...
                'Color',c_map(i_T,:),'MarkerEdgeColor', c_map(i_T,:));...,'MarkerFaceColor', c_map(i_N,:))
            set(hSCF_line(i_N,i_T), 'HandleVisibility','off');
        end
        
        
        
        str={sprintf('$T = %.2f$',T),sprintf('$\\eta = %.2f$',eta)};
        yshift = yshifts(i_T);
        xpos = xpos_T(i_T);
        ypos = ypos_T(i_T);
        xvec = [xpos, xpos];
        yvec = [ypos, ypos+yshift];
        if yshift < 0
            valign='bottom';
        else
            valign='top';
        end
        arrow_annot(i_T)=annotation('textarrow',xvec,yvec,'String',str,...
            'Units','normalized',...
            'LineWidth',.5,...
            'VerticalAlignment',valign, 'HorizontalAlignment','center',...
            'interpreter','latex',...
            'Color','black','FontSize', fontsize_annotation,...
            'HeadStyle','plain','HeadWidth',2,'HeadLength',3);
    
    end
    for i_N = 1:N_N
        dispname=sprintf('$N = (%d)^2$', sqrtN_vals(i_N));
        plot(NaN,NaN, 'DisplayName',dispname,...
            'LineStyle', 'none', 'LineWidth', 1, ...
            'Marker', symbols(i_N), 'MarkerSize', 4, ...
            'Color','k','MarkerEdgeColor', 'k');
    end





    % Legend, axes etc
    hLegend = legend('Location', 'Southwest','interpreter','latex','FontSize', 20,...
        'NumColumns',2);
    xlim([2e-2 .3]);
    ylim(ylim_vals);

    hXLabel = xlabel('$r/L$','interpreter','latex');
    hYLabel = ylabel('$L^\eta C_m(r)$','interpreter','latex');

    % Font
    set(gca, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)


    % Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'log','Yscale', 'log')
    ax_full = gca;
    if (~strcmp(label,''))
        h_text=add_subfig_label(gca,label,"se","log","log",fontsize_subfiglabels);
%         pos=get(gca,'Position');
%         dim=[pos(1) + pos(3)- labelboxdim(1) - boxbound(1),...
%             pos(2) + boxbound(2),...
%             labelboxdim];
%         annotation('textbox',dim, 'String', label, 'FitBoxToText','off',...
%             'interpreter','latex',...
%             'LineWidth', .5, ...
%             'VerticalAlignment','middle', 'HorizontalAlignment','left',...
%             'Color','black','FontSize', fontsize_annotation,...
%             'BackgroundColor','white');
    end
    figname=sprintf('%s/%s_SpinSCF_FS_3in1',basedir,curmodel);
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end
