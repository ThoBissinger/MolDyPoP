clear all
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;

% mxydata=load('mxy/rho_3.00_dynamics_LinearTime.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat'); mxyfit_FS=load('mxy/rho_3.00_DataCollapse_SCF.mat');
% xydata=load('xy/xy_dynamics_LinearTime.mat'); xyfit=load('xy/xy_CritExpFit.mat'); xyfit_FS=load('xy/xy_DataCollapse_SCF.mat');

mxydata_name='mxy/rho_3.00_dynamics_LinearTime.mat';
xydata_name='xy/xy_dynamics_LinearTime.mat';
fmxydata_name='fmxy/fmxy_dynamics_LinearTime.mat';

mxyfit_name='mxy/rho_3.00_CritExpFit.mat';
xyfit_name='xy/xy_CritExpFit.mat';
fmxyfit_name='fmxy/fmxy_CritExpFit.mat';

mxyfit_FS_name='mxy/rho_3.00_DataCollapse_SCF.mat';
xyfit_FS_name='xy/xy_DataCollapse_SCF.mat';
fmxyfit_FS_name='fmxy/fmxy_DataCollapse_SCF.mat';

mxyfit_FSMag_name='mxy/rho_3.00_etafit_FSMag.mat';
xyfit_FSMag_name='xy/xy_etafit_FSMag.mat';
fmxyfit_FSMag_name='fmxy/fmxy_etafit_FSMag.mat';
for i_model = [1]
    
    if (i_model == 1)
        curmodel="mxy";
        curtitle="";
        
        data=matfile(mxydata_name);
        fitdata=matfile(mxyfit_name);
        FS_SCF_fitdata=matfile(mxyfit_FS_name);
        FS_Mag_fitdata=matfile(mxyfit_FSMag_name);

        L_vals=[9.25,18.5,37,74,148];
        r_min = 4;
        r_max = 70;
        
        T_min_pow = 0;
        T_max_pow = .21;
        T_min_exp = .155;
        T_max_exp = .35;
        T_ticksep = .05;
%         T_select_FS=find(T_vals <= .185);
        
        % Confidence intervals for error bars
        pow_confint_val = .95;
        exp_confint_val = .95;
    elseif (i_model == 2)
        curmodel="xy";
%         curtitle="SXY model";
        curtitle="";
        
        data=matfile(xydata_name);
        fitdata=matfile(xyfit_name);
        FS_SCF_fitdata=matfile(xyfit_FS_name);
        FS_Mag_fitdata=matfile(xyfit_FSMag_name);

        
        L_vals=[16,32,64,128,256];
        r_min = 8;
        r_max = 55;
        
        
        T_min_pow = 0;
        T_max_pow = 1.8;
        T_min_exp = 1.43;
        T_max_exp = 3.00;
        T_ticksep = .5;
        
%         T_select_FS=find(T_vals <= 1.60);
        
        
        pow_confint_val = .95;
        exp_confint_val = .95;
    elseif (i_model == 3)
        curmodel="fmxy";
        curtitle="";
        
        data=matfile(fmxydata_name);
        fitdata=matfile(fmxyfit_name);
        FS_SCF_fitdata=matfile(fmxyfit_FS_name);
        FS_Mag_fitdata=matfile(fmxyfit_FSMag_name);
        
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
    end
    T_vals=data.('T_vals');
    sqrtN_vals=data.('sqrtN_vals');
    absM_av=cell2mat(data.('absM_av'));
    SCF_Spin_av = data.('SCF_Spin_av');
    rbin = data.('rbin');
    
    eta_vals = fitdata.('eta_vals');
    eta_vals_FS_SCF = FS_SCF_fitdata.('eta_vals');
    eta_vals_FS_Mag = FS_Mag_fitdata.('eta_vals');
    eta_sd = fitdata.('eta_sd');
    xi_vals = fitdata.('xi_vals');
    xi_sd = fitdata.('xi_sd');
    
    exp_fit=fitdata.('param_SCFSpin_Exp');
    pow_fit=fitdata.('param_SCFSpin_Pow');

    i_N = numel(sqrtN_vals);
    r_vals = rbin{i_N,1};
    T_select=1:numel(T_vals);
    T_select_FS=find(eta_vals_FS_SCF);


    [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% Create basic plot
    hold on
    N_sqrtN = numel(sqrtN_vals);
    N_T = numel(T_select);

    c_map = linspecer(N_T);

    figure(i_model)
    subplot(1,2,1)
%     heta_line = line(T_vals(T_select),eta_vals(i_N,T_select));
    heta_line = errorbar(T_vals(T_select),eta_vals(i_N,T_select),...
        eta_sd(i_N,T_select));
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
        xi_sd(i_N,T_select));
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


    
    figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/static/%s_EtaXi',curmodel);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    
    
    %% 2   Two in one plot
    %% 2.1 Left axis: power law fit
    fig=figure(i_model + 10);
%     set(fig,'defaultAxesColorOrder',['red'; 'blue']);
    hax=gca;
    yyaxis 'left'
    hold on;
%     heta_line = plot(T_vals(T_select),eta_vals(i_N,T_select));
%     set(heta_line, ...
%         'LineStyle', 'none', 'LineWidth', 2, ...
%         'Color','red',...
%         'Marker','*','MarkerSize',10);
    heta_line = errorbar(T_vals(T_select),eta_vals(i_N,T_select),...
        eta_sd(i_N,T_select));
    
%         eta_errors_signed(1,T_select),eta_errors_signed(2,T_select));
    set(heta_line, ...
        'LineStyle', 'none', 'LineWidth', 1, ...
        'Color','red');
    heta_FS_line = plot(T_vals(T_select_FS),eta_vals_FS_SCF(T_select_FS),'*');
    set(heta_FS_line, ...
        'LineStyle', 'none', 'LineWidth', 1, ...
        'Marker','o','MarkerSize',10,...
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
    set(gca, 'FontSize', 12)
%     set(hLegend, 'FontSize', 12)
    set([hXLabel, hYLabel_eta], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', 'red', 'YTick', 0:.05:1, ...
        'LineWidth', .5,'Xscale', 'lin','XTick',0:T_ticksep:T_max_exp)

    
%     hTBKT_line = xline(T_KT(i_N));
%     set(hTBKT_line, ...
%         'LineStyle', ':', 'LineWidth', 1.5, ...
%         'Color','green');
    hTstar_line = xline(T_star(i_N));
    set(hTstar_line, ...
        'LineStyle', ':', 'LineWidth', 2, ...
        'Color','green');
    text(T_star(i_N) - T_ticksep/6, .015,'$T^*$',...
        'VerticalAlignment','bottom','HorizontalAlignment','right',...
        'interpreter','latex',...
        'Color','green','FontSize', 12);

    %% 2.2 Right axis: exponential fit
    yyaxis 'right'
    T_select=find(T_vals > T_min_exp);
%     hxi_line = plot(T_vals(T_select),xi_vals(i_N,T_select));
%     set(hxi_line, ...
%         'LineStyle', 'none', 'LineWidth', 2, ...
%         'Color','blue',...
%         'Marker','+','MarkerSize',10);
    hxi_line = errorbar(T_vals(T_select),xi_vals(i_N,T_select),...
        xi_sd(i_N,T_select));
%         xi_errors_signed(1,T_select),xi_errors_signed(2,T_select));
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
    set(gca, 'FontSize', 12)
    
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
    
    
    figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/static/%s_EtaXi_2in1',curmodel);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
end
    
%% 3 Eta compare: Fit to (256)^2, fit to FS Magnetization, fit to FS SCF
figure(i_model + 20);
hold on
N_sqrtN = numel(sqrtN_vals);
N_T = numel(T_select);

c_map = linspecer(3);

heta_line_SCF = errorbar(T_vals,eta_vals(i_N,:),...
    eta_sd(i_N,:));
set(heta_line_SCF, ...
    'LineStyle', 'none', 'LineWidth', 1.8, ...
    'Color',c_map(1,:),...
    'DisplayName','Fit to SCF');

heta_line_FS_Mag = line(T_vals,eta_vals_FS_Mag);
set(heta_line_FS_Mag, ...
    'LineStyle', 'none', 'LineWidth', 1.5, ...
    'Color',c_map(2,:),...
    'Marker','o','MarkerSize',10,...
    'DisplayName','Fit to M');

heta_line_FS_SCF = line(T_vals(find(eta_vals_FS_SCF)),eta_vals_FS_SCF(find(eta_vals_FS_SCF)));
set(heta_line_FS_SCF, ...
    'LineStyle', 'none', 'LineWidth', 1.5, ...
    'Color',c_map(3,:),...
    'Marker','s','MarkerSize',10,...
    'DisplayName','FS Fit to SCF');

hLegend=legend('location','northwest');
xlim([T_min_pow T_max_pow]);
ylim([0 .4]);
hXLabel = xlabel('$T$','interpreter','latex');
hYLabel_eta = ylabel('$\eta$','interpreter','latex');
hTitle = title(curtitle);
set(gca, 'FontName', 'cmr12')
set([hXLabel, hYLabel_eta], 'FontName', 'cmr12')
set(gca, 'FontSize', 12)
set(hLegend, 'FontSize', 12)
set([hXLabel, hYLabel_eta], 'FontSize', 12)
set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:.05:1, ...
    'LineWidth', .5,'Xscale', 'lin','XTick',0:T_ticksep:T_max_exp)

figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/static/%s_Eta_Compare',curmodel);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end