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
        
        exp_fit=fitdata.('param_SCFSpin_Exp');
        pow_fit=fitdata.('param_SCFSpin_Pow');
        eta_vals=fitdata.('eta_vals');
        xi_vals=fitdata.('xi_vals');
        
        i_N = numel(sqrtN_vals);
        r_vals = rbin{i_N,1};
        T_select=1:numel(T_vals);
        
        L_vals=[9.25,18.5,37,74,148];
        r_min = 4;
        r_max = 70;
        
        T_min_pow = 0;
        T_max_pow = .21;
        T_min_exp = .155;
        T_max_exp = .35;
        T_ticksep = .05;
        
        % Confidence intervals for error bars
        pow_confint_val = .95;
        exp_confint_val = .95;
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
        
        exp_fit=fitdata.('param_SCFSpin_Exp');
        pow_fit=fitdata.('param_SCFSpin_Pow');
        eta_vals=fitdata.('eta_vals');
        xi_vals=fitdata.('xi_vals');
        
        i_N = numel(sqrtN_vals);
        r_vals = rbin{i_N,1};
        T_select=1:numel(T_vals);
        
        L_vals=sqrtN_vals;
        r_min = 8;
        r_max = 55;
        
        
        T_min_pow = 0;
        T_max_pow = 1.8;
        T_min_exp = 1.43;
        T_max_exp = T_vals(end);
        T_ticksep = .5;
        
        pow_confint_val = .95;
        exp_confint_val = .95;
    end


    [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% Create basic plot
    hold on
    N_sqrtN = numel(sqrtN_vals);
    N_T = numel(T_select);
    eta_vals = zeros(N_sqrtN,N_T);
    eta_errors = zeros(2,N_T);
    eta_errors_signed = zeros(2,N_T);
    xi_vals = zeros(N_sqrtN,N_T);
    xi_errors = zeros(2,N_T);
    xi_errors_signed = zeros(2,N_T);
    c_map = linspecer(N_T);
    for i = 1:N_T
        i_T = T_select(i);
        T = T_vals(i_T);
%         subplot(2,1,1)
        
        SCF=cell2mat(SCF_Spin_av(i_N,i_T));
        curpowfit=fit_PowSCF(r_vals,SCF,r_min,r_max,0);
        pow_fit{i} = curpowfit;
        eta_vals(i_N,i_T) = curpowfit.eta;
        pow_confint{i} = confint(curpowfit,.90);
        eta_errors(:,i_T) = pow_confint{i}(:,2);
        
        curexpfit=fit_ExpSCF(r_vals,SCF,r_min,r_max);
        exp_fit{i} = curexpfit;
        xi_vals(i_N,i_T) = curexpfit.xi;
        exp_confint{i} = confint(curexpfit,.90);
        xi_errors(:,i_T) = exp_confint{i}(:,2);
        
    end
    eta_errors_signed(1,:) = eta_errors(1,:) - eta_vals(i_N,:);
    eta_errors_signed(2,:) = eta_errors(2,:) - eta_vals(i_N,:);
    
    xi_errors_signed(1,:) = xi_errors(1,:) - xi_vals(i_N,:);
    xi_errors_signed(2,:) = xi_errors(2,:) - xi_vals(i_N,:);
    
    figure(i_model)
    subplot(1,2,1)
%     heta_line = line(T_vals(T_select),eta_vals(i_N,T_select));
    heta_line = errorbar(T_vals(T_select),eta_vals(i_N,T_select),...
        eta_errors_signed(1,T_select),eta_errors_signed(2,T_select));
    set(heta_line, ...
        'LineStyle', 'none', 'LineWidth', 1, ...
        'Color','red');
%         'Marker', '^', 'MarkerSize', 4, ...
%         'Color',c_map(i,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i,:))
    xlim([T_min_pow T_max_pow]);
    ylim([0 .4]);
    hXLabel = xlabel('$T$','interpreter','latex');
    hYLabel_eta = ylabel('$\eta$','interpreter','latex');
    hTitle = title(curtitle);
    set(gca, 'FontName', 'cmr12')
    set([hXLabel, hYLabel_eta], 'FontName', 'cmr12')
    set(gca, 'FontSize', 8)
%     set(hLegend, 'FontSize', 12)
    set([hXLabel, hYLabel_eta], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:.05:1, ...
        'LineWidth', .5,'Xscale', 'lin','XTick',0:T_ticksep:T_max_exp)
    
    subplot(1,2,2)
    hxi_line = errorbar(T_vals(T_select),xi_vals(i_N,T_select),...
        xi_errors_signed(1,T_select),xi_errors_signed(2,T_select));
%     line(T_vals(T_select),xi_vals(i_N,T_select));
    set(hxi_line, ...
        'LineStyle', 'none', 'LineWidth', 1, ...
        'Color','blue');
%         'Marker', 'o', 'MarkerSize', 6, ...
%         'Color',c_map(i,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i,:))
    xlim([T_min_exp T_max_exp]);
    ylim([0 .8*L_vals(i_N)]);
    hXLabel = xlabel('$T$','interpreter','latex');
    hYLabel_eta = ylabel('$\xi$','interpreter','latex');
    hTitle = title(curtitle);
    set(gca, 'FontName', 'cmr12')
    set([hXLabel, hYLabel_eta], 'FontName', 'cmr12')
    set(gca, 'FontSize', 8)
%     set(hLegend, 'FontSize', 12)
    set([hXLabel, hYLabel_eta], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'YTick', 0:floor(L_vals(i_N)/10):.8*L_vals(i_N), ...
        'LineWidth', .5,'Xscale', 'lin','XTick',0:T_ticksep:T_max_exp)


    
    figname=sprintf('plots/%s_EtaXi',curmodel);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    
    
    %% 2   Two in one plot
    %% 2.1 Left axis: power law fit
    fig=figure(i_model + 2);
%     set(fig,'defaultAxesColorOrder',['red'; 'blue']);
    hax=gca;
    yyaxis 'left'
    heta_line = errorbar(T_vals(T_select),eta_vals(i_N,T_select),...
        eta_errors_signed(1,T_select),eta_errors_signed(2,T_select));
    set(heta_line, ...
        'LineStyle', 'none', 'LineWidth', 1, ...
        'Color','red');
    hetacrit_line = yline(.25);
    set(hetacrit_line, ...
        'LineStyle', '--', 'LineWidth', 1.5, ...
        'Color','red');
    text(T_ticksep/3, .95*.25,'$\eta_C = 1/4$',...
        'VerticalAlignment','top','interpreter','latex',...
        'Color','red','FontSize', 12);
    ylim([0 .4]);
    hYLabel_eta = ylabel('$\eta$','interpreter','latex','Color','red');
    hTitle = title(curtitle);
    set(gca, 'FontName', 'cmr12')
    set([hXLabel, hYLabel_eta], 'FontName', 'cmr12')
    set(gca, 'FontSize', 8)
%     set(hLegend, 'FontSize', 12)
    set([hXLabel, hYLabel_eta], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', 'red', 'YTick', 0:.05:1, ...
        'LineWidth', .5,'Xscale', 'lin','XTick',0:T_ticksep:T_max_exp)

    
    hTBKT_line = xline(T_KT(i_N));
    set(hTBKT_line, ...
        'LineStyle', ':', 'LineWidth', 1.5, ...
        'Color','green');
    text(T_KT(i_N) - T_ticksep/6, .015,'$T_{\textrm{BKT}}$',...
        'VerticalAlignment','bottom','HorizontalAlignment','right',...
        'interpreter','latex',...
        'Color','green','FontSize', 12);

    %% 2.2 Right axis: exponential fit
    yyaxis 'right'
    hxi_line = errorbar(T_vals(T_select),xi_vals(i_N,T_select),...
        xi_errors_signed(1,T_select),xi_errors_signed(2,T_select));
    set(hxi_line, ...
        'LineStyle', 'none', 'LineWidth', 1, ...
        'Color','blue');
    hLhalf_line = yline(L_vals(i_N)/2);
    set(hLhalf_line, ...
        'LineStyle', '--', 'LineWidth', 1.5, ...
        'Color','blue');
    text(T_max_exp - T_ticksep/3, 1.05*L_vals(i_N)/2,'$L/2$',...
        'VerticalAlignment','bottom','HorizontalAlignment','right',...
        'interpreter','latex',...
        'Color','blue','FontSize', 12);
    ylim([0 .7*L_vals(i_N)]);
    hYLabel_xi = ylabel('$\xi$','interpreter','latex');
    hTitle = title(curtitle);
    set(gca, 'FontName', 'cmr12')
    set([hYLabel_xi], 'FontName', 'cmr12','Color','blue')
    set(gca, 'FontSize', 8)
    
    xlim([0 T_max_exp]);
    hXLabel = xlabel('$T$','interpreter','latex');
    set(hXLabel, 'FontName', 'cmr12')
%     set(hLegend, 'FontSize', 12)
    set([hXLabel, hYLabel_xi], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', 'blue', ...
        'YTick', 0:floor(L_vals(i_N)/10):2*L_vals(i_N), ...
        'LineWidth', .5,'Xscale', 'lin','XTick',0:T_ticksep:T_max_exp)

    
end
    