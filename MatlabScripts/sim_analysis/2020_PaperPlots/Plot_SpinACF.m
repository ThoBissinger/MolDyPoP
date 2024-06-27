clear all
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=0;

mxydata_name='mxy/rho_3.00_dynamics_LinearTime.mat';
xydata_name='xy/xy_dynamics_LinearTime.mat';
fmxydata_name='fmxy/fmxy_dynamics_LinearTime.mat';
mxyfit_name='mxy/rho_3.00_CritExpFit.mat';
xyfit_name='xy/xy_CritExpFit.mat';
fmxyfit_name='fmxy/fmxy_CritExpFit.mat';

% mxydata=load('mxy/rho_3.00_dynamics.mat'); 
% mxyfit=load('mxy/rho_3.00_CritExpFit.mat');
% xydata=load('xy/lf0_qreduced.mat'); xyfit=load('xy/lf0_CritExpFit.mat');
% /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_dynamics.mat

for i_model = 1:1

    % figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=load(mxydata_name);
        fitdata=load(mxyfit_name);
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        ACF_Spin=data.('ACF_Spin');
        TCF_times=data.('TCF_times');
        qbin = data.('qbin');
                
        i_N = numel(sqrtN_vals);
        q_vals = qbin{i_N,1};
        T_select=[9:2:13,14:21];
        T_select=1:numel(T_vals);
        
        L_vals=[9.25,18.5,37,74,148];
        r_min = 6;
        r_max = 35;
        
        T_max = .4;
        FSplot_min = .13;
        FSplot_max = .23;
        FS_Tstep = .05;
        T_offset = .0025; % For text in inset
    elseif (i_model == 2)
        curmodel="xy";
        curtitle="SXY model";
        
        data=load(xydata_name);
        fitdata=load(xyfit_name);
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
    elseif (i_model == 3)
        curmodel="fmxy";
        curtitle="FMXY model";
        
        data=load(fmxydata_name);
        fitdata=load(fmxyfit_name);
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        ACF_Spin=data.('ACF_Spin');
        TCF_times=data.('TCF_times');
        qbin = data.('qbin');
                
        i_N = numel(sqrtN_vals);
        q_vals = qbin{i_N,1};
        T_select=1:numel(T_vals);
        
        L_vals=[9.25,18.5,37,74,148];
        r_min = 6;
        r_max = 35;
        
        T_max = .4;
        FSplot_min = .13;
        FSplot_max = .23;
        FS_Tstep = .05;
        T_offset = .0025; % For text in inset
    end


%     [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% Create basic plot
    figure(i_model)
    hold on
    N_sqrtN = numel(sqrtN_vals);
    c_map = linspecer(numel(T_select));
%     coeffs = cell2mat(param_TCFSpin_omega_1_q_DO(i_N,:)');
    for i = 1:length(T_select)
%         fig1=figure(10*i_model + i);
%         hold on
        i_T = T_select(i);

        T = T_vals(i_T);
        dispname=sprintf('T = %.3f', T);
        
        ACF_Spin_cur = ACF_Spin{i_N,i_T};
        t_vals=TCF_times{i_N,i_T};
        fitfunc_DO=@(c) TCF_Spin_cur(1)*exp(-c(1) * t_vals/2) .* cos(c(2) * t_vals - c(3))/cos(c(3));
%         c_forq = coeffs_cur(3*(i_q - 1) + (1:3));
        
        plot(t_vals,ACF_Spin_cur,...
            'Color',c_map(i,:),'DisplayName',dispname,...
            'LineStyle', '--', ...
            'Marker', '^', 'MarkerSize', 4, ...
            'LineWidth',1.5);
        hold on;
%         legend show
%         legend('Location','best')
    end
    
    h_axis = gca;
    set(h_axis,'xscale','log');
        
    hXLabel = xlabel('$t$','interpreter','latex');
    hYLabel = ylabel('$C_{m}(q,t)$','interpreter','latex');




    %% Legend, axes etc
    hLegend = legend('Location', 'SouthWest','interpreter','latex','FontSize', 10,...
        'NumColumns',2);
    
    hTitle = title(curtitle);

    %% Font
    set(h_axis, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set(h_axis, 'FontSize', 12)
    set([hXLabel, hYLabel], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


    %% Adjust axes properties
    set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'log','Yscale', 'lin')
    
    figname=sprintf('plots/%s_SpinACF_sqrtN_%d',curmodel,sqrtN_vals(i_N));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
   
end
    