clear all
close all
% addpath /home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=0;

mxydata=load('mxy/rho_3.00_ACF.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat');
xydata=load('xy/lf0_eq.mat'); xyfit=load('xy/lf0_CritExpFit.mat');


    
    
    curmodel="mxy";
    curtitle="MXY model";

    data=mxydata;
    fitdata=mxyfit;
    T_vals=data.('T_vals');
    sqrtN_vals=data.('sqrtN_vals');
%     absM_av=cell2mat(data.('absM_av'));
    ACF_MSD = data.('ACF_MSD');
    ACF_anglediff = data.('ACF_anglediff');
%     rbin = data.('rbin');
    averaging_times = data.('averaging_times');
    t_vals = cell2mat(averaging_times(end,1));

    expfit=fitdata.('param_SCFSpin_Exp');

    i_N = numel(sqrtN_vals);
    i_N_plot = numel(sqrtN_vals);
    T_select=[9:2:13,14:21];
    T_max = .4;
    T_grating = 100;
    T_select=[1:find(T_vals > T_max,1) - 1];

    L_vals=[9.25,18.5,37,74,148];
    r_min = 6;
    r_max = 35;

    labels=["$N = 2^{8}$", "$N = 2^{10}$", "$N = 2^{12}$", "$N = 2^{14}$", "$N = 2^{16}$"];
%     [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% Create MSD plot
    MSDfig=figure(1);
    movegui('northwest')
    hold on
    c_map = jet(T_grating);
    for i = 1:numel(T_select)
        i_T = T_select(i);
        T = T_vals(i_T);
        color_ind = floor(T / T_max * T_grating);
        
        curMSD=cell2mat(ACF_MSD(i_N,i_T));
        hMSD_line(i_T) = line(t_vals,curMSD);
        set(hMSD_line(i_T), ...
            'LineStyle', '--', 'LineWidth', 1, ...
            'Color',c_map(color_ind,:));
        
        
        f_MSDfit = fittype('4*D*x');
        t_indices = find(t_vals < 1e3 & t_vals > 1e0);
        cur_t_vals = t_vals(t_indices);
        curMSD = curMSD(t_indices);
        c = fit(cur_t_vals(:), curMSD(:), f_MSDfit); % ,'StartPoint',c0,'Weights',weights,...
        D_r(i) = c.D;
        MSDfits{i} = c;
        MSD_D(i) = c.D;
        MSD_confint{i} = confint(c,.95);
        MSD_errors(:,i) = MSD_confint{i}(:);
        

    end
    

    %% Legend, axes etc
    MSDhcbar = colorbar;
    colormap jet
    caxis([T_vals(T_select(1)) T_vals(T_select(end))])
    
    hXLabel = xlabel('$t$','interpreter','latex');
    hYLabel = ylabel('$\langle (r(t) - r(0))^2\rangle$','interpreter','latex');
    hTitle = title(curtitle);

    %% Font
    set(gca, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set(gca, 'FontSize', 8)
    set([hXLabel, hYLabel], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


    %% Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'log','Yscale', 'log')
    MSD_ax_full = gca;

    
    %% Inset
    MSD_ax_inset = axes('Position',[.5 .15 .28 .28]);
    set(MSD_ax_inset, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'LineWidth', .5, 'Xscale', 'lin')
    hMSD_D_line = line(T_vals(T_select),D_r);
    c_map = linspecer(1);
    set(hMSD_D_line, ...
        'LineStyle', '--', 'LineWidth', 1.5, ...
        'Marker', '+', 'MarkerSize', 6, ...
        'Color',c_map(1,:))
    
%     hXLabel = xlabel('$T$','interpreter','latex');
%     hYLabel = ylabel('$D_r$','interpreter','latex');
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),'$T$',...
        'VerticalAlignment','bottom','HorizontalAlignment','right',...
        'interpreter','latex',...
        'Color','black','FontSize', 10);
    text(NW(1),NW(2),'$D_{r}$',...
        'VerticalAlignment','top','HorizontalAlignment','left',...
        'interpreter','latex',...
        'Color','black','FontSize', 10);
    set(gca, 'FontName', 'cmr12')
%     set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set([gca], 'FontSize', 8)
%     set([hXLabel, hYLabel], 'FontSize', 12)
    
    
    figname='plots/mxy_MSD';
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
   
    
    %% Create ang_MSD plot
    angMSDfig=figure(2);
    movegui('north')
    hold on
    c_map = jet(T_grating);
    for i_N = 1:numel(sqrtN_vals)
        
        for i = 1:numel(T_select)
            i_T = T_select(i);
            T = T_vals(i_T);
            color_ind = floor(T / T_max * T_grating);

            curanglediff=cell2mat(ACF_anglediff(i_N,i_T));
            if i_N == i_N_plot
                hanglediff_line(i_T) = line(t_vals,curanglediff);
                set(hanglediff_line(i_T), ...
                    'LineStyle', '--', 'LineWidth', 1, ...
                    'Color',c_map(color_ind,:));
            end

            f_angMSDfit = fittype('2*D*x');
            t_indices = find(t_vals < 1e3 & t_vals > 1e0);
            cur_t_vals = t_vals(t_indices);
            curanglediff = curanglediff(t_indices);
            c = fit(cur_t_vals(:), curanglediff(:), f_angMSDfit); % ,'StartPoint',c0,'Weights',weights,...
            D_theta(i_N,i) = c.D;
            angMSDfits{i_N,i} = c;
            angMSD_D(i_N,i) = c.D;
            angMSD_confint{i_N,i} = confint(c,.95);
            angMSD_errors(:,i_N,i) = angMSD_confint{i}(:);

        end
    end

    %% Legend, axes etc
    ang_hcbar = colorbar;
    colormap jet
    caxis([T_vals(T_select(1)) T_vals(T_select(end))])
    
    hXLabel = xlabel('$t$','interpreter','latex');
    hYLabel = ylabel('$\langle (\theta(t) - \theta(0))^2\rangle$','interpreter','latex');
    hTitle = title(curtitle);

    %% Font
    set(gca, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set(gca, 'FontSize', 8)
    set([hXLabel, hYLabel], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


    %% Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'log','Yscale', 'log')
    ang_ax_full = gca;

    
    %% Inset
    ang_ax_inset = axes('Position',[.5 .15 .28 .28]);
    set(ang_ax_inset, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'LineWidth', .5, 'Xscale', 'lin', 'Yscale', 'lin')
    c_map = linspecer(numel(sqrtN_vals));
    for i_N = 1:numel(sqrtN_vals)
        hangMSD_D_line{i_N} = line(T_vals(T_select),D_theta(i_N,:));
        
        set(hangMSD_D_line{i_N}, ...
            'LineStyle', '--', 'LineWidth', 1.5, ...
            'Marker', '+', 'MarkerSize', 6, ...
            'Color',c_map(i_N,:))
%         'DisplayName', labels(i_N), ...
    end
    ylim([0 .9]);
%     hXLabel = xlabel('$T$','interpreter','latex');
%     hYLabel = ylabel('$D_{\theta}$','interpreter','latex');
    hLegend = legend(labels(1),labels(2));
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),'$T$',...
        'VerticalAlignment','bottom','HorizontalAlignment','right',...
        'interpreter','latex',...
        'Color','black','FontSize', 10);
    text(NW(1),NW(2),'$D_{\theta}$',...
        'VerticalAlignment','top','HorizontalAlignment','left',...
        'interpreter','latex',...
        'Color','black','FontSize', 10);
    set(gca, 'FontName', 'cmr12')
%     set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set([gca], 'FontSize', 8)
%     set([hXLabel, hYLabel], 'FontSize', 12)
    hLegend = legend('Location', 'West','interpreter','latex','FontSize', 8);
    
    figname='plots/mxy_andMSD';
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
   