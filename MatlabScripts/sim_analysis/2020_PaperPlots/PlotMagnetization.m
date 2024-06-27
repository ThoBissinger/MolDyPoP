clear all
close all
% addpath /home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
saveswitch=0;

mxydata=load('mxy/rho_3.00_eq.mat');
xydata=load('xy/lf0_eq.mat');

for i_model = 1:2

    figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=mxydata;
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        absM_av=cell2mat(data.('absM_av'));
        L_vals=[9.25,18.5,37,74,148];
        T_max = .4;
        FSplot_min = .15;
        FSplot_max = .255;
        FS_Tstep = .05;
        T_offset = .0025; % For text in inset
        labels=["N = 2^{8}", "N = 2^{10}", "N = 2^{12}", "N = 2^{14}", "N = 2^{16}"];
    else
        curmodel="xy";
        curtitle="SXY model";
        
        data=xydata;
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        absM_av=cell2mat(data.('absM_av'));
        L_vals=sqrtN_vals;
        T_max = 3.5;
        FSplot_min = 1.2;
        FSplot_max = 2.1;
        FS_Tstep = .5;
        T_offset = .03; % For text in inset
        labels=['N = 2^{8}', 'N = 2^{10}', 'N = 2^{12}', 'N = 2^{14}', 'N = 2^{16}'];
    end


    [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% Create basic plot
    figure
    hold on
    N_sqrtN = numel(sqrtN_vals);
%     c_map = turbo(N_sqrtN+2); c_map= c_map(2:end-1,:);

    c_map = colormap_sqrtN();
    
    for i_N = 1 : N_sqrtN
%         hMag_line(i_N) = line(T_vals,absM_av(i_N,:));
        
%         set(hMag_line(i_N), ...
%             'LineStyle', '--', 'LineWidth', 1.5, ...
%             'Color', c_map(i_N,:))
%         set(hMag_line(i_N), ...
%             'LineStyle', 'none')
        fitfunc=mag_fitfuncs{i_N};

        TT=linspace(.7*T_star(i_N), T_C(i_N));
        hfit_line(i_N) = line(TT,fitfunc(T_C(i_N),TT));
        set(hfit_line(i_N), ...
            'LineStyle', '-', 'LineWidth', 2, ...
            'Color', c_map(i_N,:))

        
    end
    for i_N = 1 : N_sqrtN
        hMag_dots(i_N) = line(T_vals,absM_av(i_N,:));
        set(hMag_dots(i_N), ...
            'LineStyle', 'none', ...
            'Marker', 'o', 'MarkerSize', 4, ...
            'DisplayName', labels(i_N), ...
            'Color',c_map(i_N,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))
        
        hcrossover_line(i_N) = line(T_star(i_N),crossover_M(i_N));
        set(hcrossover_line(i_N), ...
            'LineStyle', 'none', ...
            'Marker', 's', 'MarkerSize', 10, ...
            'Color', c_map(i_N,:), 'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))
    end


    %% Legend, axes etc
    if (i_model == 1)
        hLegend = legend([hMag_dots(1), hMag_dots(2), hMag_dots(3), hMag_dots(4)], ...
            '$N = 2^{10}$', '$N = 2^{12}$', '$N = 2^{14}$', '$N = 2^{16}$', ...
            'Location', 'SouthWest','interpreter','latex','FontSize', 20);
    else
        hLegend = legend([hMag_dots(1), hMag_dots(2), hMag_dots(3), hMag_dots(4), hMag_dots(5)], ...
            '$N = 2^{8}$', '$N = 2^{10}$', '$N = 2^{12}$', '$N = 2^{14}$', '$N = 2^{16}$', ...
            'Location', 'SouthWest','interpreter','latex','FontSize', 20);
        
    end
    % hlegend(2).LineStyle = '-';
    % lineEntry = findobj(hLegend.EntryContainer, 'Object',hMag_dots(2));
    % entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    % entryMarker.LineWidth = .5;
    hXLabel = xlabel('$T$','interpreter','latex');
    hYLabel = ylabel('$\langle|M|\rangle$','interpreter','latex');
    hTitle = title(curtitle);
    xlim([0, T_max]);

    %% Font
    set(gca, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set(gca, 'FontSize', 8)
    set(hLegend, 'FontSize', 10)
    set([hXLabel, hYLabel], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


    %% Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:.2:1, ...
        'LineWidth', .5)
    ax_full = gca;
    
    %% Inset
    ax_inset = axes('Position',[.6 .6 .28 .28]);
    set(ax_inset, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'YTick', 0:FS_Tstep:T_max, ...
        'LineWidth', .5, 'Xscale', 'log')
    ylim([FSplot_min, FSplot_max]);
    
    hTC_line = line(sqrtN_vals.^2,T_C);
    hTstar_line = line(sqrtN_vals.^2,T_star);
    hTKT_line = line(sqrtN_vals.^2,T_KT);
    
    c_map = linspecer(3);
    set(hTC_line, ...
        'LineStyle', '--', 'LineWidth', 1.5, ...
        'Marker', '^', 'MarkerSize', 6, ...
        'Color',c_map(1,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(1,:))
    set(hTstar_line, ...
        'LineStyle', '--', 'LineWidth', 1.5, ...
        'Marker', 's', 'MarkerSize', 7, ...
        'Color',c_map(2,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(2,:))
    set(hTKT_line, ...
        'LineStyle', '--', 'LineWidth', 1.5, ...
        'Marker', 'diamond', 'MarkerSize', 6, ...
        'Color',c_map(3,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(3,:))
    
    hXLabel = xlabel('$N$','interpreter','latex');
    hYLabel = ylabel('$T$','interpreter','latex');
    set(gca, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set([gca], 'FontSize', 8)
    set([hXLabel, hYLabel], 'FontSize', 12)
    
    %% Text
    text(sqrtN_vals(2)^2,T_C(2)+T_offset,'$T_C$',...
        'VerticalAlignment','bottom','interpreter','latex',...
        'Color',c_map(1,:),'FontSize', 10);
    text(sqrtN_vals(2)^2,T_star(2)+T_offset,'$T^*$',...
        'VerticalAlignment','bottom','interpreter','latex',...
        'Color',c_map(2,:),'FontSize', 10);
    text(sqrtN_vals(2)^2,T_KT(2)-T_offset,'$T_{BKT}$',...
        'VerticalAlignment','top','interpreter','latex',...
        'Color',c_map(3,:),'FontSize', 10);
    
    figname=sprintf('plots/%s_Magnetization',curmodel);
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end
    