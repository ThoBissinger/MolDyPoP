clear all
close all
% addpath /home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=0;

mxydata=load('mxy/rho_3.00_eq.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat');
xydata=load('xy/lf0_eq.mat'); xyfit=load('xy/lf0_CritExpFit.mat');


for i_model = 1:2

    figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=mxydata;
        fitdata=mxyfit;
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        absM_av=cell2mat(data.('absM_av'));
        SCF_Spin_av = data.('SCF_Spin_av');
        rbin = data.('rbin');
        
        expfit=fitdata.('param_SCFSpin_Exp');
        powfit=fitdata.('param_SCFSpin_Pow');
        eta_vals=fitdata.('eta_vals');
        xi_vals=fitdata.('xi_vals');
        
        i_N = numel(sqrtN_vals);
        r_vals = rbin{i_N,1};
        T_select=[9:2:13,14:21];
        T_select=[9:2:23];
        
        L_vals=[9.25,18.5,37,74,148];
        r_min = 6;
        r_max = 35;
        
        T_max = .4;
        FSplot_min = .13;
        FSplot_max = .23;
        FS_Tstep = .05;
        T_offset = .0025; % For text in inset
%         labels=['N = 2^{10}', 'N = 2^{12}', 'N = 2^{14}', 'N = 2^{16}'];
    else
        curmodel="xy";
        curtitle="SXY model";
        
        data=xydata;
        fitdata=xyfit;
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        absM_av=cell2mat(data.('absM_av'));
        SCF_Spin_av = data.('SCF_Spin_av');
        rbin = data.('rbin');
        
        expfit=fitdata.('param_SCFSpin_Exp');
        powfit=fitdata.('param_SCFSpin_Pow');
        eta_vals=fitdata.('eta_vals');
        xi_vals=fitdata.('xi_vals');
        
        i_N = numel(sqrtN_vals);
        r_vals = rbin{i_N,1};
        T_select=[17:3:38];
        
        L_vals=sqrtN_vals;
        r_min = 8;
        r_max = 55;
        
        
        T_max = 3.5;
        FSplot_min = 1.2;
        FSplot_max = 2.1;
        FS_Tstep = .5;
        T_offset = .03; % For text in inset
        
    end


    [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% Create basic plot
    hold on
    N_sqrtN = numel(sqrtN_vals);
    c_map = linspecer(numel(T_select));
    for i = 1:numel(T_select)
        i_T = T_select(i);
        T = T_vals(i_T);
%         subplot(2,1,1)
        
        SCF=cell2mat(SCF_Spin_av(i_N,i_T));

        rr = linspace(r_vals(1),r_vals(end));
        curpowfit=fit_PowSCF(r_vals,SCF,r_min,r_max,0);
        powfit{i} = curpowfit;
        curexpfit=fit_ExpSCF(r_vals,SCF,r_min,r_max);
        expfit{i} = curexpfit;
        errfunc=@(r,y1,y2) sum(( r > r_min ) .* (r < r_max) .* (y1 - y2').^2);
%         if ( errfunc(r_vals,SCF,curpowfit(r_vals)) < errfunc(r_vals,SCF,curexpfit(r_vals)))
        if ( T <= T_C(i_N))
            hFit_line(i_T) = line(rr,curpowfit(rr));
            set(hFit_line(i_T), 'LineStyle', '--')
        else
            hFit_line(i_T) = line(rr,curexpfit(rr));
            set(hFit_line(i_T), 'LineStyle', ':')
        end
%         hPowFit_line(i_T) = line(rr,rr.^-eta_vals(i_N,i_T) * SCF(end)/rr(end)^-eta_vals(i_N,i_T));
        set(hFit_line(i_T), ...
            'LineWidth', 2, ...
            'HandleVisibility','Off',...
            'Color',c_map(i,:));

        dispname=sprintf('$T = %.2f$', T);
        hSCF_line(i_T) = line(r_vals(1:2:end),SCF(1:2:end));
        set(hSCF_line(i_T), ...
            'LineStyle', 'none', 'LineWidth', 1, ...
            'DisplayName',dispname,...
            'Marker', 'o', 'MarkerSize', 6, ...
            'Color',c_map(i,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i,:))
    end
    





    %% Legend, axes etc
    hLegend = legend('Location', 'Northeast','interpreter','latex','FontSize', 20,...
        'NumColumns',3);
    xlim([3 r_max]);
    ylim([1e-3 3e0]);
    legend('location','northeast');
    % hlegend(2).LineStyle = '-';
    % lineEntry = findobj(hLegend.EntryContainer, 'Object',hMag_dots(2));
    % entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    % entryMarker.LineWidth = .5;
    
    hXLabel = xlabel('$r$','interpreter','latex');
    hYLabel = ylabel('$\langle s(0)s(r)\rangle$','interpreter','latex');
    hTitle = title(curtitle);

    %% Font
    set(gca, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set(gca, 'FontSize', 8)
    set(hLegend, 'FontSize', 12)
    set([hXLabel, hYLabel], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


    %% Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'lin','Yscale', 'log')
    ax_full = gca;
    
    figname=sprintf('plots/%s_SpinSCF_sqrtN_%d',curmodel,sqrtN_vals(i_N));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
   
end
    