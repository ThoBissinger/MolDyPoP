clear all;
close all;

addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;
fig_select=[1:2];
fig_base = "/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/FFT";
models = ["mxy", "xy", "fmxy"];
% models = ["mxy", "xy"];
% curmodel = "mxy";
for i_model = 1:3
    curmodel = models(i_model);
    % mxydata=load('mxy/rho_3.00_integ.mat');
    % mxydata=load('mxy/rho_3.00_dynamics_better_q.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat');
    if (curmodel == "mxy")
        data=load('mxy/rho_3.00_dynamics_LinearTime.mat');
        fit_FS=load('mxy/rho_3.00_DataCollapse_SCF.mat');
        T_select=[2,4,8,10,11:16];
        T_select=[1:17];
    %     T_select=[17:19];
        xmax = .45;
        L_vals = [9.25,18.5,37,74,148];
        z=1;
        n_q = 1;
    elseif (curmodel == "xy")
        data=load('xy/xy_dynamics_LinearTime.mat');
        fit_FS=load('xy/xy_DataCollapse_SCF.mat');
        T_select=[8,10,13,16,18:24,26];
        xmax=2.5 ;
        L_vals=[16,32,64,128,256];
        z=1;
        n_q = 1;
    elseif (curmodel == "fmxy")
        data=load('fmxy/fmxy_dynamics_LinearTime.mat');
        fit_FS=load('fmxy/fmxy_DataCollapse_SCF.mat');
%         T_select=[2,4,8,10,11:16];
        T_select=[1:17];
        xmax = .45;
        L_vals = [9.25,18.5,37,74,148];
        z=1;
        n_q = 1;
    end

%     plot_type={"w"};
% 
%     gxx=data.('gxx');
%     gyy=data.('gyy');
%     gmparmpar=data.('gmparmpar');
    gmperpmperp=data.('gmperpmperp');
%     gtt=data.('gtt');
%     gww=data.('gww');
%     chimxq_av=data.('chimxq_av');
%     chimyq_av=data.('chimyq_av');
%     chimparq_av=data.chimparq_av;
%     chimperpq_av=data.chimperpq_av;
    % TODO chiteq. Missing in current dataset
    qbin=data.('qbin');
    averaging_times=data.('averaging_times');

    eta_vals=fit_FS.('eta_vals');
    sqrtN_vals=data.('sqrtN_vals');
    T_vals=data.('T_vals');

    %% CALCULATION
    if (curmodel == "mxy")
        c_vel_mxy = zeros(size(T_select));
    elseif (curmodel == "xy")    
        c_vel_xy = zeros(size(T_select));
    elseif  (curmodel == "fmxy")
        c_vel_fmxy = zeros(size(T_select));
    end
    i_N=5;
    q_vals=qbin{i_N,1};
    om_mperp = zeros(numel(T_select),numel(q_vals));
    c_map = linspecer(numel(T_select));
    for ind_T = 1:numel(T_select)
        sqrtN=sqrtN_vals(i_N);
        i_T = T_select(ind_T);
        T=T_vals(i_T);
        dispname = sprintf('T = %.3f', T);
%             eta=eta_vals(i_T);
        
        q_select = 1:numel(q_vals);
        
        for i_q = q_select
            q = q_vals(i_q);
            q_indices = (i_q):length(q_vals):numel(gmperpmperp{i_N,i_T});
            [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                gmperpmperp{i_N,i_T}(q_indices), 0);
            [fmax,i_om_max] = max(ft_vals);
            om_mperp(ind_T,i_q) = abs(om_vals(i_om_max));
        end
        
        fit_depth = 10;
        if (curmodel == "mxy")
            P = polyfit(q_vals(1:fit_depth),om_mperp(ind_T,1:fit_depth),1);
            c_vel_mxy(ind_T) = P(1);
        elseif (curmodel == "xy")
            P = polyfit(q_vals(1:fit_depth),om_mperp(ind_T,1:fit_depth),1);
            c_vel_xy(ind_T) = P(1);
        elseif  (curmodel == "fmxy")
            P = polyfit(q_vals(1:fit_depth),om_mperp(ind_T,1:fit_depth),1);
            c_vel_fmxy(ind_T) = P(1);
        end
%         avg_depth = 5;
%         if (curmodel == "mxy")
% %             First order: (om_mperp(ind_T,2) - om_mperp(ind_T,1))/(q_vals(2) - q_vals(1));
% %             Also 1st:    (om_mperp(ind_T,4) - om_mperp(ind_T,3))/(q_vals(4) - q_vals(3));
% %                Possibly better results?
% %             Second order: (4 * om_mperp(ind_T,2) - 3 * om_mperp(ind_T,1) - om_mperp(ind_T,3))/(q_vals(2) - q_vals(1)) * .5;
% %                Doesn't work, since not equidistant
%             om_diffs_mxy(ind_T,:) = om_mperp(ind_T,2:end) - om_mperp(ind_T,1:end-1);
%             q_diffs_mxy = q_vals(2:end) - q_vals(1:end-1);
%             c_vel_mxy(ind_T) = mean(om_diffs_mxy(ind_T,1:avg_depth)./q_diffs_mxy(1:avg_depth));
%         elseif (curmodel == "xy")
%             om_diffs_xy(ind_T,:) = om_mperp(ind_T,2:end) - om_mperp(ind_T,1:end-1);
%             q_diffs_xy = q_vals(2:end) - q_vals(1:end-1);
%             c_vel_xy(ind_T) = mean(om_diffs_xy(ind_T,1:avg_depth)./q_diffs_xy(1:avg_depth));
%         elseif  (curmodel == "fmxy")
%             om_diffs_fmxy(ind_T,:) = om_mperp(ind_T,2:end) - om_mperp(ind_T,1:end-1);
%             q_diffs_fmxy = q_vals(2:end) - q_vals(1:end-1);
%             c_vel_fmxy(ind_T) = mean(om_diffs_fmxy(ind_T,1:avg_depth)./q_diffs_fmxy(1:avg_depth));
%         end
    end
    %% FIGURE 1: Dispersion
    
    if(ismember(1,fig_select))
        figure(1 + 1e4*i_model);
        i_N=5;
        q_vals=qbin{i_N,1};
        c_map = linspecer(numel(T_select));
        for ind_T = 1:numel(T_select)
            sqrtN=sqrtN_vals(i_N);
            i_T = T_select(ind_T);
            T=T_vals(i_T);
            dispname = sprintf('T = %.3f', T);
%             eta=eta_vals(i_T);
            
            plot(q_vals,om_mperp(ind_T,:),...
                'Color',c_map(ind_T,:),'DisplayName',dispname, ...
                'LineStyle', 'none', 'LineWidth',1.5, ...
                'Marker', '.', 'MarkerSize', 12);
            hold on;
    end
    
%             xlim([0,xmax]);
        h_axis = gca;
        hXLabel = xlabel('$q$','interpreter','latex');
        hYLabel = ylabel('$\omega_{m\perp}(q)$','interpreter','latex');



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Legend, axes etc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hLegend = legend('Location', 'Southeast','interpreter','latex','FontSize', 12,...
            'NumColumns',2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Font
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(h_axis, 'FontName', 'cmr12')
        set([hXLabel, hYLabel], 'FontName', 'cmr12')
        set(h_axis, 'FontSize', 14)
        set([hXLabel, hYLabel], 'FontSize', 14)
    %     set(hTitle, 'FontSize', 14, 'FontWeight' , 'bold')


        annotation_str = {sprintf('$N = (%d)^2$',sqrtN)};
        dim=[.16 .8 .1 .2];
        annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
            'interpreter','latex',...
            'VerticalAlignment','middle', 'HorizontalAlignment','left',...
            'Color','Black','FontSize', 12,...
            'BackgroundColor','white');

        figname=sprintf('%s/Dispersion/%s_omega_0_mperp_sqrtN_%d',fig_base,curmodel,sqrtN);
        fprintf('Creating figure %s\n',figname)
        if(saveswitch == 1)
           exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
    end

    
end
%% FIGURE 2: SW velocity
if(ismember(2,fig_select))
    figure(2);
    i_N=5;
    h_axis = gca;
    plot(T_vals(T_select),c_vel_mxy,...
        'Color','red','DisplayName','c (MXY)', ...
        'LineStyle', '-', 'LineWidth',1.5, ...
        'Marker', '.', 'MarkerSize', 12);
    hold on;
    plot(T_vals(T_select),c_vel_fmxy,...
        'Color','blue','DisplayName','c (FMXY)', ...
        'LineStyle', '-', 'LineWidth',1.5, ...
        'Marker', '.', 'MarkerSize', 12);
    hXLabel = xlabel('$T$','interpreter','latex');
    hYLabel = ylabel('$c_{m\perp}$','interpreter','latex');



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Legend, axes etc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hLegend = legend('Location', 'Northeast','interpreter','latex','FontSize', 12,...
        'NumColumns',2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Font
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set(h_axis, 'FontSize', 14)
    set([hXLabel, hYLabel], 'FontSize', 14)
%     set(hTitle, 'FontSize', 14, 'FontWeight' , 'bold')


    annotation_str = {sprintf('$N = (%d)^2$',sqrtN)};
    dim=[.16 .1 .1 .2];
    annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
        'interpreter','latex',...
        'VerticalAlignment','middle', 'HorizontalAlignment','left',...
        'Color','Red','FontSize', 14,...
        'BackgroundColor','white');

    figname=sprintf('%s/Dispersion/sound_velcotiy_mperp_fmxy_mxy_sqrtN_%d',fig_base,sqrtN);
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end

