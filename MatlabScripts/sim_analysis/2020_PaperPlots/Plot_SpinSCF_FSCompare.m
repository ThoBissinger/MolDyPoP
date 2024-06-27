clear all
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=0;

mxydata=load('mxy/rho_3.00_eq.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat'); mxyfit_FS=load('mxy/rho_3.00_DataCollapse_SCF.mat');
% xydata=load('xy/lf0_eq.mat');
xydata=load('xy/xy_dynamics_LinearTime.mat');
xyfit=load('xy/lf0_CritExpFit.mat');
xyfit_FS=load('xy/xy_DataCollapse_SCF.mat');

pathbase="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots";
for i_model = 2:2

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
        eta_vals_FS = mxyfit_FS.('eta_vals');
        xi_vals=fitdata.('xi_vals');
        
        i_N = numel(sqrtN_vals);
%         r_vals = rbin{i_N,1};
        T_select=[1,3,6,11,15,16,17,18,19,21,23,26];
        T_select=[1:19];
        
        L_vals=[9.25,18.5,37,74,148];
%         r_min = 6;
%         r_max = 60;
        
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
        eta_vals_FS = xyfit_FS.('eta_vals');
        xi_vals=fitdata.('xi_vals');
        
        i_N = numel(sqrtN_vals);
%         r_vals = rbin{i_N,1};
        T_select=[3,7,11,12,13,17,19,22,25,28,30,33,35];
        T_select=1:30;
        
        L_vals=sqrtN_vals;
%         r_min = 8;
%         r_max = 55;
        
        
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
    c_map = linspecer(numel(sqrtN_vals));
    for i = 1:numel(T_select)
        i_T = T_select(i);
        T = T_vals(i_T);
        eta=eta_vals_FS(i_T);
        figure(i_model*100 + i)
        for i_N = 1:numel(sqrtN_vals)
            sqrtN = sqrtN_vals(i_N);
            SCF=abs(cell2mat(SCF_Spin_av(i_N,i_T)));
            
            r_vals = rbin{i_N,1};
            r_min = r_vals(2);
            r_max = r_vals(end);

            rr = linspace(r_vals(1),r_vals(end));
%             curpowfit=fit_PowSCF(r_vals,SCF,r_min,r_max,0);
%             powfit{i} = curpowfit;
%                 hFit_line(i_N) = line(rr,curpowfit(rr));
%                 set(hFit_line(i_N), 'LineStyle', '--')
%             set(hFit_line(i_N), ...
%                 'LineWidth', 2, ...
%                 'HandleVisibility','Off',...
%                 'Color',c_map(i_N,:));

            dispname=sprintf('$N = %d^2$', sqrtN);
            hSCF_line(i_N) = line(r_vals(1:end)/L_vals(i_N),SCF(1:end)*L_vals(i_N)^eta);
            set(hSCF_line(i_N), ...
                'LineStyle', 'none', 'LineWidth', 1, ...
                'DisplayName',dispname,...
                'Marker', '.', 'MarkerSize', 10, ...
                'Color',c_map(i_N,:),'MarkerEdgeColor', c_map(i_N,:),'MarkerFaceColor', c_map(i_N,:))
        end
    





        %% Legend, axes etc
        hLegend = legend('Location', 'Southwest','interpreter','latex','FontSize', 20,...
            'NumColumns',1);
        xlim([0 .3]);
%         ylim([1e-2 3e0]);
        % hlegend(2).LineStyle = '-';
        % lineEntry = findobj(hLegend.EntryContainer, 'Object',hMag_dots(2));
        % entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
        % entryMarker.LineWidth = .5;

        hXLabel = xlabel('$r/L$','interpreter','latex');
        hYLabel = ylabel('$L^\eta C_m(r)$','interpreter','latex');
%         hTitle = title(curtitle);

        %% Font
        set(gca, 'FontName', 'cmr12')
        set([hXLabel, hYLabel], 'FontName', 'cmr12')
        set(gca, 'FontSize', 8)
        set(hLegend, 'FontSize', 12)
        set([hXLabel, hYLabel], 'FontSize', 12)
%         set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


        %% Adjust axes properties
        set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
            'LineWidth', .5,'Xscale', 'log','Yscale', 'lin')
        ax_full = gca;
        
        text(.05,.5,sprintf('$T = %.3f$',T),...
            'Units','normalized',...
            'VerticalAlignment','top', 'HorizontalAlignment','left',...
            'interpreter','latex',...
            'Color','Red','FontSize', 10);
        text(.05,.45,sprintf('$\\eta = %.2f$',eta),...
            'Units','normalized',...
            'VerticalAlignment','top', 'HorizontalAlignment','left',...
            'interpreter','latex',...
            'Color','Red','FontSize', 10);

        figname=sprintf('%s/plots/%s_SpinSCF_FS_T_%.3f',pathbase,curmodel,T);
        fprintf('Creating figure %s\n',figname)
        if(saveswitch == 1)
           exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
    
    end
end
    