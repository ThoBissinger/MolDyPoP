clear all;
close all;

addpath(genpath('/home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;
fig_select=[1:3];
fig_base = "/home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/FFT";
models = ["mxy", "xy", "fmxy"];
% models = ["mxy", "xy"];
% curmodel = "mxy";
for i_model = 1:1
    curmodel = models(i_model);
    % mxydata=load('mxy/rho_3.00_integ.mat');
    % mxydata=load('mxy/rho_3.00_dynamics_better_q.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat');
    if (curmodel == "mxy")
        data=load('mxy/rho_3.00_dynamics_LinearTime.mat');
        fit_FS=load('mxy/rho_3.00_DataCollapse_SCF.mat');
        T_select=[2,4,8,10,11:19];
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
        T_select=[2,4,8,10,11:19];
        xmax = .45;
        L_vals = [9.25,18.5,37,74,148];
        z=1;
        n_q = 1;
    end

    plot_type={"w", "mperp", "te","mpar", "m"};
    plot_type={"mperp", "te","mpar", "m"};
%     plot_type={"w"};

    gxx=data.('gxx');
    gyy=data.('gyy');
    gmparmpar=data.('gmparmpar');
    gmperpmperp=data.('gmperpmperp');
    gtt=data.('gtt');
    gww=data.('gww');
    chimxq_av=data.('chimxq_av');
    chimyq_av=data.('chimyq_av');
    chimparq_av=data.chimparq_av;
    chimperpq_av=data.chimperpq_av;
    % TODO chiteq. Missing in current dataset
    qbin=data.('qbin');
    averaging_times=data.('averaging_times');

    eta_vals=fit_FS.('eta_vals');
    sqrtN_vals=data.('sqrtN_vals');
    T_vals=data.('T_vals');

    for i_plot_type = 1:numel(plot_type)
        plot_type_cur = plot_type{i_plot_type};
        
        
        if(ismember(1,fig_select))
    
            for i_T = T_select
                figure(i_T + 1e4*i_model + 1e5 * i_plot_type);
                i_N=5;
    
                sqrtN=sqrtN_vals(i_N);
                T=T_vals(i_T);
                eta=eta_vals(i_T);
    
    
                q_vals=qbin{i_N,i_T};
                q_select = 1:6;
                q_offset = 0;
                c_map = linspecer(numel(q_select));
                for i_q = q_select
                    q = q_vals(i_q + q_offset);
                    dispname=sprintf('$q = %.3f$', q);
                    q_indices = (i_q+q_offset):length(q_vals):numel(gxx{i_N,i_T});
                    if (plot_type_cur == "m")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gxx{i_N,i_T}(q_indices) + gyy{i_N,i_T}(q_indices), 0);
                    elseif (plot_type_cur == "mpar")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gmparmpar{i_N,i_T}(q_indices), 0);
                    elseif (plot_type_cur == "mperp")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gmperpmperp{i_N,i_T}(q_indices), 0);
                    elseif (plot_type_cur == "te")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gtt{i_N,i_T}(q_indices), 0);
                    elseif (plot_type_cur == "w")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gww{i_N,i_T}(q_indices), 0);
                    end
                    plot(om_vals/q,abs(ft_vals)*q^(3-eta),...
                        'Color',c_map(i_q,:),'DisplayName',dispname, ...
                        'LineStyle', '-', 'LineWidth',1.5, ...
                        'Marker', 'none', 'MarkerSize', 5);
                    hold on;
                end
                xlim([0,xmax]);
                h_axis = gca;
                hXLabel = xlabel('$\omega / q$','interpreter','latex');
                if (plot_type_cur == "m")
                    hYLabel = ylabel('$S_{m}(q,\omega) q^{3-\eta}$','interpreter','latex');
                elseif (plot_type_cur == "mpar")
                    hYLabel = ylabel('$S_{m\parallel}(q,\omega) q^{3-\eta}$','interpreter','latex');
                elseif (plot_type_cur == "mperp")
                    hYLabel = ylabel('$S_{m\perp}(q,\omega) q^{3-\eta}$','interpreter','latex');
                elseif (plot_type_cur == "te")
                    hYLabel = ylabel('$S_{\theta}(q,\omega) q^{3-\eta}$','interpreter','latex');
                elseif (plot_type_cur == "w")
                    hYLabel = ylabel('$S_{w}(q,\omega) q^{3-\eta}$','interpreter','latex');
                end
                
    
    
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Legend, axes etc
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hLegend = legend('Location', 'Northeast','interpreter','latex','FontSize', 14,...
                'NumColumns',1);
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Font
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                set(h_axis, 'FontName', 'cmr12')
                set([hXLabel, hYLabel], 'FontName', 'cmr12')
                set(h_axis, 'FontSize', 14)
                set([hXLabel, hYLabel], 'FontSize', 14)
            %     set(hTitle, 'FontSize', 14, 'FontWeight' , 'bold')
    
    
                annotation_str = {sprintf('$N = (%d)^2$',sqrtN),sprintf('$T = %.3f$',T),sprintf('$\\eta = %.2f$',eta)};
                dim=[.16 .7 .1 .2];
                annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                    'interpreter','latex',...
                    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
                    'Color','Red','FontSize', 14,...
                    'BackgroundColor','white');
    
                figname=sprintf('%s/S_by_om_q_scaling/%s_S_%s_sqrtN_%d_T_%.3f',fig_base,curmodel,plot_type_cur,sqrtN,T);
                fprintf('Creating figure %s\n',figname)
                if(saveswitch == 1)
                   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
                end
            end
        end
    
    
    
    
        if(ismember(2,fig_select))
    
            i_N_start = 3;
            q = qbin{i_N_start,1}(1);
            for i_T = T_select
                figure(i_T + 100 + 1e4*i_model + 1e5 * i_plot_type);
                T=T_vals(i_T);
                eta=eta_vals(i_T);
                c_map = linspecer(numel(i_N_start : 5));
                for i_N = i_N_start:5
                    sqrtN=sqrtN_vals(i_N);
                    L=L_vals(i_N);
    
                    q_vals=qbin{i_N,i_T};
                    i_q = find(min(abs(q_vals - q)) == abs(q_vals - q));
    
                    dispname=sprintf('$N = (%d)^2$', sqrtN);
                    q_indices = (i_q):length(q_vals):numel(gxx{i_N,i_T});
                    if (plot_type_cur == "m")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gxx{i_N,i_T}(q_indices) + gyy{i_N,i_T}(q_indices), 0);
                    elseif (plot_type_cur == "mpar")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gmparmpar{i_N,i_T}(q_indices), 0);
                    elseif (plot_type_cur == "mperp")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gmperpmperp{i_N,i_T}(q_indices), 0);
                    elseif (plot_type_cur == "te")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gtt{i_N,i_T}(q_indices), 0);
                    elseif (plot_type_cur == "w")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gww{i_N,i_T}(q_indices), 0);
                    end
                    % Rough estimate of the integral
                    S_q = sum(abs(ft_vals(1:end))) * (om_vals(2) - om_vals(1));
        %             S_q = chimxq_av{i_N,i_T}(i_q) + chimyq_av{i_N,i_T}(i_q);
    
                    plot(om_vals,abs(ft_vals) / S_q,...
                        'Color',c_map(i_N - i_N_start + 1,:),'DisplayName',dispname, ...
                        'LineStyle', '-', 'LineWidth',1.5, ...
                        'Marker', 'none', 'MarkerSize', 5);
                    hold on;
                end
                xlim([0,xmax * q]);
                h_axis = gca;
                hXLabel = xlabel('$\omega$','interpreter','latex');
                if (plot_type_cur == "m")
                    hYLabel = ylabel('$S_{m}(q,\omega) / S_{m}(q)$','interpreter','latex');
                elseif (plot_type_cur == "mpar")
                    hYLabel = ylabel('$S_{m\parallel}(q,\omega) / S_{m\parallel}(q)$','interpreter','latex');
                elseif (plot_type_cur == "mperp")
                    hYLabel = ylabel('$S_{m\perp}(q,\omega) / S_{m\perp}(q)$','interpreter','latex');
                elseif (plot_type_cur == "te")
                    hYLabel = ylabel('$S_{\theta}(q,\omega) / S_{\theta}(q)$','interpreter','latex');
                elseif (plot_type_cur == "w")
                    hYLabel = ylabel('$S_{w}(q,\omega) / S_{w}(q)$','interpreter','latex');
                end
                
    
    
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Legend, axes etc
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hLegend = legend('Location', 'Northeast','interpreter','latex','FontSize', 14,...
                'NumColumns',1);
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Font
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                set(h_axis, 'FontName', 'cmr12')
                set([hXLabel, hYLabel], 'FontName', 'cmr12')
                set(h_axis, 'FontSize', 14)
                set([hXLabel, hYLabel], 'FontSize', 14)
            %     set(hTitle, 'FontSize', 14, 'FontWeight' , 'bold')
    
    
                annotation_str = {sprintf('$T = %.3f$',T),sprintf('$q = %.3f$',q),...
                    sprintf('$\\eta = %.2f$',eta)};
                dim=[.70 .49 .1 .2];
                annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                    'interpreter','latex',...
                    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
                    'Color','Red','FontSize', 14,...
                    'BackgroundColor','white');
    
                figname=sprintf('%s/S_vs_om_fixed_q/%s_S_%s_FS_T_%.3f',fig_base,curmodel,plot_type_cur,T);
                fprintf('Creating figure %s\n',figname)
                if(saveswitch == 1)
                   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
                end
            end
        end
    
    
    
        if(ismember(3,fig_select))
    
            i_N_start = 1;
            q = qbin{i_N_start,1}(1);
            for i_T = T_select
                figure(i_T + 300 + 1e4*i_model + 1e5 * i_plot_type);
                T=T_vals(i_T);
                eta=eta_vals(i_T);
                c_map = linspecer(numel(i_N_start : 5));
                for i_N = i_N_start:5
                    sqrtN=sqrtN_vals(i_N);
                    L=L_vals(i_N);
    
                    q_vals=qbin{i_N,i_T};
                    i_q = n_q;
    
                    dispname=sprintf('$N = (%d)^2$', sqrtN);
                    q_indices = (i_q):length(q_vals):numel(gxx{i_N,i_T});
                    if (plot_type_cur == "m")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gxx{i_N,i_T}(q_indices) + gyy{i_N,i_T}(q_indices), 0);
                    elseif (plot_type_cur == "mpar")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gmparmpar{i_N,i_T}(q_indices), 0);
                    elseif (plot_type_cur == "mperp")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gmperpmperp{i_N,i_T}(q_indices), 0);
                    elseif (plot_type_cur == "te")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gtt{i_N,i_T}(q_indices), 0);
                    elseif (plot_type_cur == "w")
                        [ft_vals,om_vals]=FT_correlation(averaging_times{i_N,i_T},...
                            gww{i_N,i_T}(q_indices), 0);
                    end
                    % Rough estimate of the integral
                    S_q = sum(abs(ft_vals(1:end))) * (om_vals(2) - om_vals(1));
        %             S_q = chimxq_av{i_N,i_T}(i_q) + chimyq_av{i_N,i_T}(i_q);
    
                    plot(om_vals*L^z,abs(ft_vals)/ L^z / S_q,...
                        'Color',c_map(i_N - i_N_start + 1,:),'DisplayName',dispname, ...
                        'LineStyle', '-', 'LineWidth',1.5, ...
                        'Marker', 'none', 'MarkerSize', 5);
                    hold on;
                end
                xlim([0,xmax * q * L_vals(i_N_start)^z]);
                h_axis = gca;
                hXLabel = xlabel('$\omega L^z$','interpreter','latex');
                if (plot_type_cur == "m")
                    hYLabel = ylabel('$S_{m}(q,\omega) / L^z S_{m}(q)$','interpreter','latex');
                elseif (plot_type_cur == "mpar")
                    hYLabel = ylabel('$S_{m\parallel}(q,\omega) / L^z S_{m\parallel}(q)$','interpreter','latex');
                elseif (plot_type_cur == "mperp")
                    hYLabel = ylabel('$S_{m\perp}(q,\omega) / L^z S_{m\perp}(q)$','interpreter','latex');
                elseif (plot_type_cur == "te")
                    hYLabel = ylabel('$S_{\theta}(q,\omega) / L^z S_{\theta}(q)$','interpreter','latex');
                elseif (plot_type_cur == "w")
                    hYLabel = ylabel('$S_{w}(q,\omega) / L^z S_{w}(q)$','interpreter','latex');
                end
                
    
    
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Legend, axes etc
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hLegend = legend('Location', 'Northeast','interpreter','latex','FontSize', 14,...
                'NumColumns',1);
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Font
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                set(h_axis, 'FontName', 'cmr12')
                set([hXLabel, hYLabel], 'FontName', 'cmr12')
                set(h_axis, 'FontSize', 14)
                set([hXLabel, hYLabel], 'FontSize', 14)
            %     set(hTitle, 'FontSize', 14, 'FontWeight' , 'bold')
    
    
                annotation_str = {sprintf('$T = %.3f$',T),sprintf('$n_q = %d$',i_q),...
                    sprintf('$\\eta = %.2f$',eta),sprintf('$z = %.2f$',z)};
                dim=[.70 .37 .1 .2];
                annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                    'interpreter','latex',...
                    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
                    'Color','Red','FontSize', 14,...
                    'BackgroundColor','white');
    
                figname=sprintf('%s/SbyL_vsomL/%s_S_%s_FS_n_q_%d_z_%.2f_T_%.3f',fig_base,curmodel,plot_type_cur,n_q,z,T);
                fprintf('Creating figure %s\n',figname)
                if(saveswitch == 1)
                   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
                end
            end
        end
    end
end