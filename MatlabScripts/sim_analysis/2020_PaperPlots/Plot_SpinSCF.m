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
label='$(a)$';
labelboxdim=[.05 .06];
boxbound=[.01 .01];
for i_model = 1

    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        sqrtN=256;
        data_dirs={"211201_LongerTime/mxy_3.00", "210715_LinearTimeSampling/mxy_3.00"};
        sampfilenames={"samp_Dynamics", "samp_Dynamics"};
        T_vals=[.16 .17 .18 .19 .20 .21];
        T_dirs=["T_.16" "T_.17" "T_.18" "T_.19" "T_.20" "T_.21"];
        r_min = 2;
        r_max = 60;
       
    else
        curmodel="xy";
        curtitle="SXY model";
        
        r_min = 3;
        r_max = 55;
    end
    N_T = numel(T_vals);

    %% 1 Gather data
    data_rbin=cell(1,N_T);
    data_cf=cell(1,N_T);
    for i_T = 1:N_T
        T_dir=T_dirs(i_T);
        i_dir=1;
        curfile=sprintf('%s/%s/sqrtN_%d/%s/%s.mat',data_rootdir,data_dirs{i_dir},sqrtN,T_dir,sampfilenames{i_dir});
        while (~ isfile(curfile))
            i_dir = i_dir + 1;
            curfile=sprintf('%s/%s/sqrtN_%d/%s/%s.mat',data_rootdir,data_dirs{i_dir},sqrtN,T_dir,sampfilenames{i_dir});
        end
        load(curfile,'rbin','SCF_Spin_av');
        data_rbin{i_T} = rbin;
        data_cf{i_T} = SCF_Spin_av;
    end

    %% Create plot
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    hold on

    c_map = turbo(N_T+1); c_map=c_map(2:end,:);
    for i_T = 1:N_T
        T = T_vals(i_T);

        r_vals = data_rbin{i_T};
        SCF= data_cf{i_T};
        rr = linspace(r_vals(1),r_vals(end));
        curpowfit=fit_PowSCF(r_vals,SCF,r_min,r_max,0);
        powfit{i_T} = curpowfit;
        curexpfit=fit_ExpSCF(r_vals,SCF,r_min,r_max);
        expfit{i_T} = curexpfit;
        errfunc=@(r,y1,y2) sum(( r > r_min ) .* (r < r_max) .* log((y1 - y2').^2));
        if ( errfunc(r_vals,SCF,curpowfit(r_vals)) < errfunc(r_vals,SCF,curexpfit(r_vals)))
%         if ( T <= T_C(i_N))
            hFit_line(i_T) = line(rr,curpowfit(rr));
            set(hFit_line(i_T), 'LineStyle', '-.')
        else
            hFit_line(i_T) = line(rr,curexpfit(rr));
            set(hFit_line(i_T), 'LineStyle', ':')
        end
%         hPowFit_line(i_T) = line(rr,rr.^-eta_vals(i_N,i_T) * SCF(end)/rr(end)^-eta_vals(i_N,i_T));
        set(hFit_line(i_T), ...
            'LineWidth', 2, ...
            'HandleVisibility','Off',...
            'Color',c_map(i_T,:));

        dispname=sprintf('$T = %.2f$', T);
        hSCF_line(i_T) = line(r_vals(1:2:end),SCF(1:2:end));
        set(hSCF_line(i_T), ...
            'LineStyle', 'none', 'LineWidth', .4, ...
            'DisplayName',dispname,...
            'Marker', 'o', 'MarkerSize', 4, ...
            'Color',c_map(i_T,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_T,:))
    end
    





    % Legend, axes etc
    hLegend = legend('Location', 'Northeast','interpreter','latex','FontSize', 20,...
        'NumColumns',2);
    xlim([1 r_max]);
    ylim([1e-2 3e0]);
    legend('location','northeast');
    % hlegend(2).LineStyle = '-';
    % lineEntry = findobj(hLegend.EntryContainer, 'Object',hMag_dots(2));
    % entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    % entryMarker.LineWidth = .5;
    
    hXLabel = xlabel('$r$','interpreter','latex');
    hYLabel = ylabel('$C_m(r)$','interpreter','latex');
%     hTitle = title(curtitle);

    % Font
    set(gca, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)


    % Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'lin','Yscale', 'log')
    ax_full = gca;


    % Annotations
    if (~strcmp(label,''))
        h_text=add_subfig_label(gca,label,"se","lin","log",fontsize_subfiglabels);

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
    
    if(saveswitch == 1)
        figname=sprintf('%s/%s_SpinSCF_sqrtN_%d',basedir,curmodel,sqrtN);
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
   
end
    