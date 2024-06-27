clear all
close all
% addpath /home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=0;

mxydata=load('mxy/rho_3.00_qreduced.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat');
xydata=load('xy/lf0_qreduced.mat'); xyfit=load('xy/lf0_CritExpFit.mat');


for i_model = 1:2

    figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=mxydata;
        fitdata=mxyfit;
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        gxx=data.('gxx');
        TCF_times=data.('TCF_times');
        qbin = data.('qbin');
                
        i_N = numel(sqrtN_vals);
        q_vals = qbin{i_N,1};
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


%     [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% Create basic plot
    hold on
    N_sqrtN = numel(sqrtN_vals);
    c_map = linspecer(numel(T_select));
    i_q = 3;
    c_map = linspecer(4);
%     coeffs = cell2mat(param_TCFSpin_q_DO(i_N,:)');
    for i = 1:length(T_select)
        i_T = T_select(i);
        q = q_vals(i_q);
        dispname=sprintf('q = %.3f', q);
%         coeffs_cur=param_TCFSpin_q_DO{i_N,i_T};
        T = T_vals(i_T);
        gxx_cur = gxx{i_N,i_T};
        t_vals=TCF_times{i_N,i_T};
        TCF_Spin_cur = gxx_cur(i_q:length(q_vals):end);
        fitfunc_DO=@(c) TCF_Spin_cur(1)*exp(-c(1) * t_vals/2) .* cos(c(2) * t_vals - c(3))/cos(c(3));
%         c_forq = coeffs_cur(3*(i_q - 1) + (1:3));
        
        fig1=figure(10*i_model + i);
        h = axes;
                
        plot(t_vals,real(TCF_Spin_cur),...
            'Color',c_map(1,:),'DisplayName','Simulated',...
            'LineStyle', 'none', ...
            'Marker', 'x', 'MarkerSize', 10, ...
            'LineWidth',1.5);
        hold on;
%         plot(t_vals,fitfunc_DO(c_forq),'-',...
%             'Color',c_map(2,:),'DisplayName','Fitted',...
%             'LineWidth',1.3);
        hold off;
        set(h,'xscale','log');
        
%         title(sprintf('Damped oscillator fit to time correlation function.\n T = %.2f, q = %.3f, %s', T, q, modeldata));
        
        hXLabel = xlabel('$t$','interpreter','latex');
        hYLabel = ylabel('$C_{m}(q,t)$','interpreter','latex');
%         legend show
%         legend('Location','best')
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
    
    figname=sprintf('plots/%s_SpinTCF_sqrtN_%d',curmodel,sqrtN_vals(i_N));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
   
end
    