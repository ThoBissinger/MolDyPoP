run initialization_script;
basedir=sprintf('%s/plots/TCF/ExamplesWithFit',fig_base);
saveswitch=1;

for i_model = [1,3,4]

    if (i_model == 1)
        curmodel="mxy";
        curmodelcaps="MXY";
        curtitle="MXY model";
        runmax=500;
        
        data=matfile(mxydata_AdjustedTime_name);
        FS_Mag_fitdata=matfile(mxyfit_FSMag_name);
        fitdata_FS=matfile(mxyfit_FS_name);
        fitdata_TCF=matfile(mxyfit_TCF_q_name);

        T_select = [4,5,6,7];
        
        L_vals=[9.25,18.5,37,74,148];
        t_max_vals = [5e2, 8e2, 2e3, 2e3, 6e3];
        markerspacing = [20,20,20,20,13];

        t_max_vals_short = [5e2, 8e2, 2e3, 5e2, 1e3];
        markerspacing_short = [8,8,8,8,4];

        om_max_max=.045;

        sqrtN_select=[4,5];
    elseif (i_model == 2)
        curmodel="xy";
        curmodelcaps="XY";
        curtitle="SXY model";
        runmax=250;
        
        data=matfile(xydata_name);
        FS_Mag_fitdata=matfile(xyfit_FSMag_name);
        fitdata_FS=matfile(xyfit_FS_name);
        fitdata_TCF=matfile(xyfit_TCF_q_name);

        om_max_max=.045;

        L_vals=[16 32 64 128 256];
    elseif (i_model == 3)
        curmodel="fmxy";
        curmodelcaps="FMXY";
        curtitle="FMXY model";
        runmax=500;
        
        data=matfile(fmxydata_AdjustedTime_name);
        FS_Mag_fitdata=matfile(fmxyfit_FSMag_name);
        fitdata_FS=matfile(fmxyfit_FS_name);
        fitdata_TCF=matfile(fmxyfit_TCF_q_name);

        T_select = [4,5,6,7];
        
        L_vals=[9.25,18.5,37,74,148];
        t_max_vals = [1e3, 2e3, 3e3, 5e3, 1e4];
        markerspacing = [5,5,5,5,20];
        
        t_max_vals_short = [5e2, 8e2, 2e3, 5e2, 1e3];
        markerspacing_short = [4,4,4,4,4];
        
        om_max_max=.045;

        sqrtN_select=[4,5];
    elseif (i_model == 4)
        curmodel="xy_s";
        curmodelcaps="XY";
        curtitle="SXY model";
        runmax=250;
        
        data=matfile(xysdata_AdjustedTime_name);
        FS_Mag_fitdata=matfile(xysfit_FSMag_name);
        fitdata_FS=matfile(xysfit_FSMag_name);
        fitdata_TCF=matfile(xysfit_TCF_q_name);

        T_select = [1,3,5];

        L_vals=[16 32 64 128];
        t_max_vals = [2e3, 4e3, 1e4, 1.5e4, 1e4];
        markerspacing = [5,5,5,5];
        
        t_max_vals_short = [5e2, 8e2, 2e3, 5e2, 1e3];
        markerspacing_short = [1,1,1,1,1];
        
        om_max_max=.045;

        sqrtN_select=[3,4];
    end
    res_factor=4;
    res_function=@(t,tau) resolution_Laplace_pleateau(t,tau,tau);
    q_select = [1,3,6];
    marker_types=["v" "^" "o" "d"];
    absM_av=cell2mat(data.('absM_av'));
    absM_var=cell2mat(data.('absM_var'));
    eta_vals_FS_Mag = FS_Mag_fitdata.('eta_vals');

%     plot_type=["m", "mperp", "mpar", "w"];
%     plot_type=["m", "mperp", "mpar", "w" "te"];
    plot_type=["mperp"];
    if (ismember("m",plot_type))
        gxx=data.('gxx');
        gyy=data.('gyy');
        chimxq_av=data.('chimxq_av');
        chimyq_av=data.('chimyq_av');
    end
    if (ismember("mpar",plot_type))
        gmparmpar=data.('gmparmpar');
        chimparq_av=data.chimparq_av;
    end
    if (ismember("mperp",plot_type))
        gmperpmperp=data.('gmperpmperp');
        chimperpq_av=data.chimperpq_av;
    end
    if (ismember("te",plot_type))
        gtt=data.('gtt');
        % TODO chiteq. Missing in current dataset
    end
    if (ismember("w",plot_type))
        gww=data.('gww');
        chiwq_av=data.('chiwq_av');
    end
    
    qbin=data.('qbin');
    averaging_times=data.('averaging_times');

    eta_vals=fitdata_FS.('eta_vals');
    sqrtN_vals=data.('sqrtN_vals');
    T_vals=data.('T_vals');
    N_sqrtN=numel(sqrtN_vals);
    N_T=numel(T_select);

    param_TCFSpin_omega_1_q_DO = fitdata_TCF.('param_TCFSpin_omega_1_q_DO');
    c_vals = fitdata_TCF.('omega_1_a');

    

    [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% 1 Create basic plot
%     i_N = 5;
%     sqrtN = sqrtN_vals(i_N);
    
    c_map = linspecer(N_T);
    for i_plot = 1:numel(plot_type)
        plot_type_cur=plot_type(i_plot);
        for ind_q = 1:numel(q_select)
            i_q = q_select(ind_q);
            for i_N = sqrtN_select
                sqrtN = sqrtN_vals(i_N);
                t_max = t_max_vals(i_N);
    
                figure;
                hold on
                
                for i = 1 : N_T
                    i_T=T_select(i);
                    T=T_vals(i_T);
                    t=averaging_times{i_N,i_T};
                    n_t=numel(t);
                    q_vals=qbin{i_N,i_T};
                    n_q=numel(q_vals);
    %                 i_q = 1;
                    q = q_vals(i_q);
                
                    dispname=sprintf('$T = %.3f$', T);
                    q_indices = (i_q):n_q:n_q*n_t;
                    
                    gamma_cur = param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+1);
                    omega_1_cur = param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+2);
                    omega_0_cur = sqrt(omega_1_cur^2 + .25 * gamma_cur^2);
                    
                    if (plot_type_cur == "m")
                        cf = gxx{i_N,i_T}(q_indices) + gyy{i_N,i_T}(q_indices);
                    elseif (plot_type_cur == "mpar")
                        cf = gmparmpar{i_N,i_T}(q_indices);
                    elseif (plot_type_cur == "mperp")
                        cf = gmperpmperp{i_N,i_T}(q_indices);
                    elseif (plot_type_cur == "te")
                        cf = gtt{i_N,i_T}(q_indices);
                    elseif (plot_type_cur == "w")
                        cf = gww{i_N,i_T}(q_indices);
                    end
                    cf= real(cf)/real(cf(1));
                    offset= -2.5*(i - 1);
                    c = fit_DampedOscillator_RealSpace(t,cf,n_period,weightexp,'omega_1');
                    
                    yline(offset,'--',...
                        'HandleVisibility','off');
    
                    h_fit_data(i) = plot(t,fitfunc_DO(t,1,c)+offset,...
                        'LineStyle', '-', ...
                        'LineWidth', 2, ...
                        'HandleVisibility', 'off', ...
                        'Color',c_map(i,:));
            
                    h_sim_data(i) = plot(t,cf + offset,...
                        'LineStyle', 'none', ...
                        'Marker', marker_types(i), 'MarkerSize', 2, ...
                        'HandleVisibility', 'off', ...
                        'Color',c_map(i,:),...
                        'MarkerFaceColor',c_map(i,:),'MarkerEdgeColor','k', ...
                        'MarkerIndices',1:markerspacing(i_N):numel(t));
            
                    h_leg_data(i) = plot(NaN*t,NaN*t,...
                        'LineStyle', '-', ...
                        'LineWidth', 1, ...
                        'Marker', marker_types(i), 'MarkerSize', 4, ...
                        'DisplayName', dispname, ...
                        'MarkerFaceColor',c_map(i,:),'MarkerEdgeColor','k', ...
                        'Color',c_map(i,:));
                    
                end
                ylim([offset-1.2 8.6])
                xlim([0 t_max]);
                hLegend = legend('Location', 'NorthEast','interpreter','latex',...
                    'NumColumns',1);
                
                hXLabel = xlabel('$t$','interpreter','latex');
                if (plot_type_cur == "m")
                    hYLabel = ylabel('$C_{m}(q,t) / \chi_{m}(q)$','interpreter','latex');
                elseif (plot_type_cur == "mpar")
                    hYLabel = ylabel('$C_{m\parallel}(q,t) / \chi_{m\parallel}(q)$','interpreter','latex');
                elseif (plot_type_cur == "mperp")
                    hYLabel = ylabel('$C_{m\perp}(q,t) / \chi_{m\perp}(q)$','interpreter','latex');
                elseif (plot_type_cur == "te")
                    hYLabel = ylabel('$C_{\theta}(q,t) / \chi_{\theta}(q)$','interpreter','latex');
                elseif (plot_type_cur == "w")
                    hYLabel = ylabel('$C_{w}(q,t) / \chi_{w}(q)$','interpreter','latex');
                end
            
                % Font
                set_fonts_default;
            
                % Adjust axes properties
                set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
                    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
                    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', [], ...
                    'LineWidth', .5)
                ax_full = gca;
                
                set(gcf,'units','centimeters','OuterPosition',[0 0 2/3*columnwidth_cm .8*columnwidth_cm]);
    
                annotation_str = {curmodelcaps,sprintf('$q = %.3f$',q)};
    %             innerpos=get(gcf,'InnerPosition');
                dim=[.16, .71, .1 .2];
                annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                    'interpreter','latex',...
                    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
                    'Color','black','FontSize', fontsize_annotation,...
                    'BackgroundColor','white');
            
                figname=sprintf('%s/%s/%s_FitExample_sqrtN_%d_q_%.3f',basedir,plot_type_cur,curmodel,sqrtN,q);
                if(saveswitch == 1)
                    fprintf('Creating figure %s\n',figname)
                    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
                    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
                end
        
                for i = 1 : N_T
                    i_T=T_select(i);
                    t=averaging_times{i_N,i_T};
                    n_t=numel(t);
                    set(h_sim_data(i),'MarkerIndices',1:markerspacing_short(i_N):n_t);
                end
                xlim([0 t_max_vals_short(i_N)]);
                figname=sprintf('%s/%s/%s_FitExample_ShortTime_sqrtN_%d_q_%.3f',basedir,plot_type_cur,curmodel,sqrtN,q);
                if(saveswitch == 1)
                    fprintf('Creating figure %s\n',figname)
                    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
                    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
                end
    
                %% Plot of FFT
                figure;
                hold on;
                om_max=0;
                gamma_max=0;
                for i = 1 : N_T
                    i_T=T_select(i);
                    T=T_vals(i_T);
                    t=averaging_times{i_N,i_T};
                    n_t=numel(t);
                    q_vals=qbin{i_N,i_T};
                    n_q=numel(q_vals);
                    q = q_vals(i_q);

                    
                
                    dispname=sprintf('$T = %.3f$', T);
                    q_indices = (i_q):n_q:n_q*n_t;
                    
                    gamma_cur = param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+1);
                    omega_1_cur = param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+2);
                    omega_0_cur = sqrt(omega_1_cur^2 + .25 * gamma_cur^2);
                    
                    if (plot_type_cur == "m")
                        cf = gxx{i_N,i_T}(q_indices) + gyy{i_N,i_T}(q_indices);
                    elseif (plot_type_cur == "mpar")
                        cf = gmparmpar{i_N,i_T}(q_indices);
                    elseif (plot_type_cur == "mperp")
                        cf = gmperpmperp{i_N,i_T}(q_indices);
                    elseif (plot_type_cur == "te")
                        cf = gtt{i_N,i_T}(q_indices);
                    elseif (plot_type_cur == "w")
                        cf = gww{i_N,i_T}(q_indices);
                    end
                    cf= real(cf)/real(cf(1));
                    c = fit_DampedOscillator_RealSpace(t,cf,n_period,weightexp,'omega_1');
                    gamma = c(1);
                    omega_0 = sqrt(c(2)^2 - gamma^2/4);
                    
                    tau=res_factor/gamma_cur;
                    res_vals=res_function(t,tau);
%                     res_vals = resolution_Quartic(t,5/gamma);
                    [ft_vals,om_vals]=FT_correlation(t, cf .* res_vals,0);
                    ft_vals = real(ft_vals);
    
                    c_fft = fit_DampedOscillator_FourierSpace(om_vals,ft_vals,cf(1));
                    om_max = max(om_max,c_fft.omega_0);
                    gamma_max = max(gamma_max,c_fft.gamma);
    %                 c = fit_DampedOscillator_RealSpace(t,cf,n_period,weightexp,'omega_9');
    
                    
                    
                    h_ft_fit_data(i) = plot(om_vals,fitfunc_DO_reciprocal(om_vals,cf(1),[c_fft.gamma,c_fft.omega_0]),...
                        'LineStyle', '-', ...
                        'LineWidth', 2, ...
                        'HandleVisibility', 'off', ...
                        'Color',c_map(i,:));
            
                    h_ft_sim_data(i) = plot(om_vals,ft_vals,...
                        'LineStyle', 'none', ...
                        'Marker', marker_types(i), 'MarkerSize', 2, ...
                        'HandleVisibility', 'off', ...
                        'Color',c_map(i,:),...
                        'LineWidth',.5,...
                        'MarkerFaceColor',c_map(i,:),'MarkerEdgeColor','k', ...
                        'MarkerIndices',1:3:numel(om_vals));
    
    %                 h_ft_sim_data(i) = plot(om_vals,ft_vals,...
    %                     'LineStyle', '--', ...
    %                     'HandleVisibility', 'off', ...
    %                     'Color',c_map(i,:));
            
                    h_leg_data(i) = plot(NaN*t,NaN*t,...
                        'LineStyle', '-', ...
                        'LineWidth', 1, ...
                        'Marker', marker_types(i), 'MarkerSize', 4, ...
                        'DisplayName', dispname, ...
                        'MarkerFaceColor',c_map(i,:),'MarkerEdgeColor','k', ...
                        'Color',c_map(i,:));
                    
                end
                old_ylim = get(gca,'ylim');
                ylim([0, 1.5*old_ylim(2)]);
                xlim([0 2.5*om_max]);
                hLegend = legend('Location', 'NorthEast','interpreter','latex',...
                    'NumColumns',1);
                
                hXLabel = xlabel('$\omega$','interpreter','latex');
                if (plot_type_cur == "m")
                    hYLabel = ylabel('$S_{m}(q,\omega) / \chi_{m}(q)$','interpreter','latex');
                elseif (plot_type_cur == "mpar")
                    hYLabel = ylabel('$S_{m\parallel}(q,\omega) / \chi_{m\parallel}(q)$','interpreter','latex');
                elseif (plot_type_cur == "mperp")
                    hYLabel = ylabel('$S_{m\perp}(q,\omega) / \chi_{m\perp}(q)$','interpreter','latex');
                elseif (plot_type_cur == "te")
                    hYLabel = ylabel('$S_{\theta}(q,\omega) / \chi_{\theta}(q)$','interpreter','latex');
                elseif (plot_type_cur == "w")
                    hYLabel = ylabel('$S_{w}(q,\omega) / \chi_{w}(q)$','interpreter','latex');
                end
            
                % Font
                set_fonts_default;
            
                % Adjust axes properties
                set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
                    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
                    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
                    'LineWidth', .5,...
                    'XScale','lin','YScale','lin')
                ax_full = gca;
                
                set(gcf,'units','centimeters','OuterPosition',[0 0 2/3*columnwidth_cm .8*columnwidth_cm]);
    
                annotation_str = {curmodelcaps,sprintf('$q = %.3f$',q)};
    %             innerpos=get(gcf,'InnerPosition');
                dim=[.21, .71, .1 .2];
                annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                    'interpreter','latex',...
                    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
                    'Color','black','FontSize', fontsize_annotation,...
                    'BackgroundColor','white');
            
                figname=sprintf('%s/%s/%s_FitExample_FFT_sqrtN_%d_q_%.3f',basedir,plot_type_cur,curmodel,sqrtN,q);
                if(saveswitch == 1)
                    fprintf('Creating figure %s\n',figname)
                    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
                    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
                end
        
                set(gca,'YScale','log');
                ylim([1e-1,300*old_ylim(2)])
                xlim([0 2.5*om_max]);
                figname=sprintf('%s/%s/%s_FitExample_FFT_logscale_sqrtN_%d_q_%.3f',basedir,plot_type_cur,curmodel,sqrtN,q);
                if(saveswitch == 1)
                    fprintf('Creating figure %s\n',figname)
                    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
                    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
                end
            end
        end
    end
end