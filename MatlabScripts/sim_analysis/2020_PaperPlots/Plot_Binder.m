run initialization_script;
basedir=sprintf('%s/plots/static',fig_base);
saveswitch=1;

for i_model = [4]

    figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=load(mxydata_name);
        L_vals=[9.25,18.5,37,74,148];
        
        
    elseif (i_model == 2)
        curmodel="xy";
        curtitle="SXY model";
        
        data=load(xydata_name);
        L_vals=[16 32 64 128 256];
    elseif (i_model == 3)
        curmodel="fmxy";
        curtitle="FMXY model";
        
        data=load(fmxydata_name);
        L_vals=[9.25,18.5,37,74,148];
    elseif (i_model == 4)
        curmodel="xy_s";
        curtitle="XY S model";
        
        data=load(xysdata_name);
        L_vals=[16 32 64 128 256];
    end
    T_vals=data.('T_vals');
    sqrtN_vals=data.('sqrtN_vals');
    Binder_cum=cell2mat(data.('Binder_cum'));
    
    labels=['N = (16)^{2}', 'N = (32)^{2}', 'N = (64)^{2}', 'N = (128)^{2}', 'N = (256)^{2}'];

    %% Create basic plot
    figure(i_model)
    hold on
    N_sqrtN = numel(sqrtN_vals);
    c_map = linspecer(N_sqrtN);
    for i_N = 1 : N_sqrtN

%         hBinder_line(i_N) = line(T_vals,Binder_cum(i_N,:));
%         set(hBinder_line(i_N), ...
%             'LineStyle', '--', 'LineWidth',2,...
%             'Marker', 'o', 'MarkerSize', 1, ...
%             'HandleVisibility', 'off', ...
%             'Color',c_map(i_N,:));

        hBinder(i_N) = line(T_vals,Binder_cum(i_N,:));
        set(hBinder(i_N), ...
            'LineStyle', '--', 'LineWidth',2,...
            'Marker', 'o', 'MarkerSize', 6, ...
            'DisplayName', labels(i_N), ...
            'Color',c_map(i_N,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))

%         hBinder_dots(i_N) = line(T_vals,Binder_cum(i_N,:));
%         set(hBinder_dots(i_N), ...
%             'LineStyle', 'none', ...
%             'Marker', 'o', 'MarkerSize', 6, ...
%             'HandleVisibility', 'off', ...
%             'Color',c_map(i_N,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))
%         
%         hcrossover_line(i_N) = line(T_star(i_N),crossover_M(i_N));
%         set(hcrossover_line(i_N), ...
%             'LineStyle', 'none', ...
%             'Marker', 's', 'MarkerSize', 10, ...
%             'Color', c_map(i_N,:), 'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))
    end


    %% Legend, axes etc
        hLegend = legend([hBinder(1), hBinder(2), hBinder(3), hBinder(4), hBinder(5)], ...
            '$N = (16)^{2}$', '$N = (32)^{2}$', '$N = (64)^{2}$', '$N = (128)^{2}$', '$N = (256)^{2}$', ...
            'Location', 'SouthWest','interpreter','latex','FontSize', 20);
    % hlegend(2).LineStyle = '-';
    % lineEntry = findobj(hLegend.EntryContainer, 'Object',hMag_dots(2));
    % entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    % entryMarker.LineWidth = .5;
    hXLabel = xlabel('$T$','interpreter','latex');
    hYLabel = ylabel('$U_{\textrm{Binder}}$','interpreter','latex');
    hTitle = title(curtitle);
%     xlim([0, T_max]);

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
    
    
    figname=sprintf('%s/%s_BinderCum',basedir,curmodel);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end
    