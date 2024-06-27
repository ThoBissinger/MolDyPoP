%% Initialization
clear
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/TCF/ShortTimeResolution',fig_base);

section_select = [3];

for i_model = [1]
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";

        data=matfile(mxydata_AdjustedTime_name);
        
        fitdata=matfile(mxyfit_name);
        fitdata_FS=matfile(mxyfit_FS_name);
        fitdata_TCF=matfile(mxyfit_TCF_q_name);

        lintime_dataset_name='dynamics_mxy_LinearTime';
        
       
%         T_select=[1:7];
%         T_select=[1,2,4,5,7];
        T_select=[1 2 6 9 11];
        L_vals=[9.25,18.5,37,74,148];
        
        y_offsets_fit=[.97,.9,.87];

        t_max_long = 8;
        t_max_short = 1;

        om_by_cq_max = 5;
        
    elseif (i_model == 2)
        curmodel="xy";
        curtitle="XY model";

        data=matfile(xydata_LepriRuffo_extended_name);

        fitdata=matfile(xyfit_name);
        fitdata_FS=matfile(xyfit_FS_name);
        fitdata_TCF=matfile(xyfit_TCF_q_name);
                
        lintime_dataset_name='dynamics_xy_LinearTime';

        T_select=[8:16];
        L_vals=[16, 32, 64, 128, 256];

    elseif (i_model == 3)
        curmodel="fmxy";
        curtitle="FMXY model";

        data=matfile(fmxydata_AdjustedTime_name);

        fitdata=matfile(fmxyfit_name);
        fitdata_FS=matfile(fmxyfit_FS_name);
        fitdata_TCF=matfile(fmxyfit_TCF_q_name);
        
        lintime_dataset_name='dynamics_fmxy_LinearTime';

        T_select=[1,2,4,5,7];
        L_vals=[9.25,18.5,37,74,148];

%         res_factor=.5;
%         res_function = resolution_Gauss;
        
    elseif (i_model == 4)
        curmodel="xy_s";
        curtitle="XY S model";

        data=matfile(xysdata_AdjustedTime_name);

        fitdata=matfile(xysfit_FSMag_name);
        fitdata_FS=matfile(xysfit_FSMag_name);
        fitdata_TCF=matfile(xysfit_TCF_q_name);
        
        lintime_dataset_name='dynamics_xy_s_AdjustedTime';

        T_select=1:5;
        L_vals=[16, 32, 64, 128];
        
    end
    res_factor=30;
    res_function = resolution_Gauss;
%     res_function=@(t,tau) resolution_Laplace_pleateau(t,tau,tau);

    i_N_start = 2;
    i_N_max = 4;
    n_trunc = 4; % For NF Psi
%     plot_type=["w", "mperp", "te","mpar", "m"];
%     plot_type=["mperp", "te","mpar", "m"];
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
    if (ismember("t",plot_type))
        gtt=data.('gtt');
        % TODO chiteq. Missing in current dataset
    end
    if (ismember("w",plot_type))
        gww=data.('gww');
        chiwq_av=data.('chiwq_av');
    end
    
    qbin=data.('qbin');
    averaging_times=data.('averaging_times');
    absM_av=data.('absM_av');


    eta_vals=fitdata_FS.('eta_vals');
    sqrtN_vals=data.('sqrtN_vals');
    T_vals=data.('T_vals');

    param_TCFSpin_omega_1_q_DO = fitdata_TCF.('param_TCFSpin_omega_1_q_DO');
    c_vals = fitdata_TCF.('omega_1_a');
    
%     res_factor=10;
    z=1;
    
    gui_locations={'northwest', 'north', 'northeast','west','center','east','southwest','south','southeast'};
    
    for i_plot_type = 1:numel(plot_type)
        plot_type_cur = plot_type(i_plot_type);
        %% Section 1 Plot to study the Nelson-Fisher high-frequency behavior
%         c_map = linspecer(numel(sqrtN_vals));
        if (ismember(1,section_select))
            c_map = linspecer(i_N_max - i_N_start +1);
%             n_q_vals=1;
%             q = qbin{i_N_start,1}(n_q_vals(1));
            q = qbin{i_N_start,1}(1);
            for i_T = T_select
                figure;
                T=T_vals(i_T);
                i_T_lintime=find_T_index(lintime_dataset_name,T);
                absM_vec=cell2mat(absM_av(:,i_T));
                eta_fitob=fit_eta_Magnetization_FS(absM_vec,L_vals);
                eta = eta_fitob.eta;
                
                for i_N = i_N_start:i_N_max
                    sqrtN=sqrtN_vals(i_N);
                    L=L_vals(i_N);
    
                    t=averaging_times{i_N,i_T};
    
    
                    q_vals=qbin{i_N,i_T};
                    i_q = find(q_vals == q);
%                     i_q = n_q_vals(i_N - i_N_start + 1);
    
                    dispname=sprintf('$N = (%d)^2$', sqrtN);
                    q_indices = (i_q):length(q_vals):length(q_vals)*length(t);
%                     gamma_cur = param_TCFSpin_omega_1_q_DO{i_N,i_T_lintime}(2*i_q+1);
%                     omega_1_cur = param_TCFSpin_omega_1_q_DO{i_N,i_T_lintime}(2*i_q+2);
%                     omega_0_cur = sqrt(omega_1_cur^2 + .25 * gamma_cur^2);
                    
                    
    
%                     t=averaging_times{i_N,i_T};

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
                    c = fit_DampedOscillator_RealSpace(t,cf,n_period,weightexp,'omega_1');
                    gamma_cur = c(1);
                    omega_1_cur = c(2);
                    omega_0_cur = sqrt(c(2)^2 - gamma_cur^2/4);

                    tau=res_factor/gamma_cur;
                    resolution_vals=resolution_Laplace(t,tau);
                    
%                     cf=fitfunc_DO(t,1,[1e-3,1e-2]);
                    [ft_vals,om_vals]=FT_correlation(t, cf .* resolution_vals, 1e5);
                    % Rough estimate of the integral
                    ft_vals=real(ft_vals);
%                     S_q = sum(real(ft_vals(1:end))) * (om_vals(2) - om_vals(1));
                    S_q = cf(1);
                    [ft_max,i_ft_max] = max(ft_vals);
        %             S_q = chimxq_av{i_N,i_T}(i_q) + chimyq_av{i_N,i_T}(i_q);
    
                    loglog(om_vals / om_vals(i_ft_max), real(ft_vals)/ L^z / S_q,...
                        'Color',c_map(i_N - i_N_start + 1,:),'DisplayName',dispname, ...
                        'LineStyle', 'none', 'LineWidth',2, ...
                        'Marker', '.', 'MarkerSize', 5);
                    hold on;
                end
%                 xlim([.1*omega_0_cur*L^z,1e2*omega_0_cur*L^z]);
                xlim_vec=[.3 6];
                y_vals=logspace(log10(xlim_vec(1)),log10(xlim_vec(2)),1e2);
                NF_Psi = NelsonFisher_Psi(y_vals,eta,n_trunc);
%                 loglog(om_vals,2*pi^2*eta^2*om_vals.^(-3-eta),...
                loglog(y_vals,NF_Psi,...
                        'Color','Black','DisplayName','NF', ...
                        'LineStyle', '--', 'LineWidth',2, ...
                        'Marker', 'none', 'MarkerSize', 5);

                xlim(xlim_vec);
                ylim([.3*min(NF_Psi),10*max(NF_Psi)]);
                h_axis = gca;
                hXLabel = xlabel('$\omega / cq$','interpreter','latex');
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
                hLegend = legend('Location', 'Northeast','interpreter','latex',...
                    'NumColumns',1);
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Font
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                set_fonts_default
    
                set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
                    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
                    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'XTick', 1:10, ...
                    'LineWidth', .5)

                annotation_str = {sprintf('$T = %.3f$',T),sprintf('$q = %.3f$',q),... sprintf('$n_q = %d$',i_q)
                    sprintf('$\\eta = %.2f$',eta),sprintf('$z = %.2f$',z)};
                dim=[.15 .15 .1 .2];
                annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                    'interpreter','latex',...
                    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
                    'Color','black','FontSize', fontsize_annotation,...
                    'BackgroundColor','white');
    
                set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

%                 figname=sprintf('%s/HighFreq/%s/%s_S_FS_n_q_%d_z_%.2f_T_%.3f_Nstart_%d_resfac_%d',basedir,plot_type_cur,curmodel,n_q,z,T,i_N_start,res_factor);
                figname=sprintf('%s/HighFreq/%s/%s_S_FS_q_%.3f_T_%.3f_resfac_%d',basedir,plot_type_cur,curmodel,q,T,res_factor);
    %                 figname=sprintf('%s/SbyL_vsomL/%s_S_%s_FS_n_q_%d_z_%.2f_tau_%de3_T_%.3f',basedir,curmodel,plot_type_cur,n_q,z,tau_laplace/1e3,T);
                fprintf('Creating figure %s\n',figname)
                if(saveswitch == 1)
                   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
                end
            end
        end

        %% Section 2 Plots of combined TCFs (log time scale and other)
        c_map = linspecer(numel(sqrtN_vals));
        if (ismember(2,section_select))
            i_N_start = 1;
            i_N_end = 2;
            n_q=1;
            q = qbin{i_N_start,1}(n_q);
            for i_T = T_select
                figure;
                T=T_vals(i_T);
                i_T_lintime=find_T_index(lintime_dataset_name,T);
                eta=eta_vals(i_T_lintime);
                c_map = linspecer(numel(i_N_start : 5));
                for i_N = i_N_start:i_N_end             
                    sqrtN=sqrtN_vals(i_N);
                    L=L_vals(i_N);
    
                    tlog=averaging_times{i_N,i_T};
                    tlin=averaging_times_lintime{i_N,i_T_lintime};
    
    
                    q_vals=qbin{i_N,i_T};
                    q_vals_lintime=qbin_lintime{i_N,i_T_lintime};
                    i_q = n_q;
    
                    dispname=sprintf('$N = (%d)^2$', sqrtN);
                    q_indices = (i_q):length(q_vals):length(q_vals)*length(tlog);
                    q_indices_lintime = (i_q):length(q_vals_lintime):length(q_vals_lintime)*length(tlin);
%                     gamma_cur = param_TCFSpin_omega_1_q_DO{i_N,i_T_lintime}(2*i_q+1);
%                     omega_1_cur = param_TCFSpin_omega_1_q_DO{i_N,i_T_lintime}(2*i_q+2);
%                     omega_0_cur = sqrt(omega_1_cur^2 + .25 * gamma_cur^2);
                    
                    tau=res_factor/gamma_cur;
    
%                     t=averaging_times{i_N,i_T};

                    if (plot_type_cur == "m")
                        cf = gxx{i_N,i_T}(q_indices) + gyy{i_N,i_T}(q_indices);
                        cf_lintime = gxx_lintime{i_N,i_T_lintime}(q_indices_lintime) ...
                            + gyy_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "mpar")
                        cf = gmparmpar{i_N,i_T}(q_indices);
                        cf_lintime = gmparmpar_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "mperp")
                        cf = gmperpmperp{i_N,i_T}(q_indices);
                        cf_lintime = gmperpmperp_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "te")
                        cf = gtt{i_N,i_T}(q_indices);
                        cf_lintime = gtt_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "w")
                        cf = gww{i_N,i_T}(q_indices);
                        cf_lintime = gww_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    end
                    [t,cf] = combine_cf(tlog,cf,tlin,cf_lintime);
%                     resolution_vals=resolution_Laplace(t,tau);

%                     [ft_vals,om_vals]=FT_correlation(t, cf .* resolution_vals, 1e6);
                    % Rough estimate of the integral
%                     S_q = sum(real(ft_vals(1:end))) * (om_vals(2) - om_vals(1));
        %             S_q = chimxq_av{i_N,i_T}(i_q) + chimyq_av{i_N,i_T}(i_q);
    
                    plot(t/L^z,real(cf / cf(1)),...
                        'Color',c_map(i_N - i_N_start + 1,:),'DisplayName',dispname, ...
                        'LineStyle', '-', 'LineWidth',2, ...
                        'Marker', 'none', 'MarkerSize', 5);
                    hold on;
                end
                xlim([0,min(max(t),4*tau)/L^z]);
                h_axis = gca;
                hXLabel = xlabel('$t / L^z$','interpreter','latex');
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
                
    
    
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Legend, axes etc
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hLegend = legend('Location', 'Northeast','interpreter','latex','FontSize', fontsize_annotation,...
                    'NumColumns',1);
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Font
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)
    
                annotation_str = {sprintf('$T = %.3f$',T),sprintf('$n_q = %d$',i_q),...
                    sprintf('$\\eta = %.2f$',eta),sprintf('$z = %.2f$',z)};
                dim=[.685 .18 .1 .2];
                annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                    'interpreter','latex',...
                    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
                    'Color','black','FontSize', fontsize_annotation,...
                    'BackgroundColor','white');
    
                set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

                figname=sprintf('%s/LogAndLin/%s/%s_S_SmallN_n_q_%d_T_%.3f',basedir,plot_type_cur,curmodel,n_q,T);
    %                 figname=sprintf('%s/SbyL_vsomL/%s_S_%s_FS_n_q_%d_z_%.2f_tau_%de3_T_%.3f',basedir,curmodel,plot_type_cur,n_q,z,tau_laplace/1e3,T);
                fprintf('Creating figure %s\n',figname)
                if(saveswitch == 1)
                   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
                end
            end
        end

        %% Section 3 Plot to study the Nelson-Fisher high-frequency with c, like Evertz-Landau fig 13
%         c_map = linspecer(numel(sqrtN_vals)+2);
        c_map = linspecer(numel(i_N_start : 5)+2);
        if (ismember(3,section_select))
            i_N_start = 3;
            q = qbin{i_N_start,1}(1);
            for i_T = T_select
                figure;
                T=T_vals(i_T);
                i_T_lintime=find_T_index(lintime_dataset_name,T);
%                 eta=eta_vals(i_T_lintime);
                M_vec = cell2mat(absM_av(:,i_T));
                eta_fitob = fit_eta_Magnetization_FS(M_vec,L_vals);
                eta = eta_fitob.eta;
                c_cur = c_vals(i_T_lintime);
                
                for i_N = i_N_start:numel(sqrtN_vals)
                    sqrtN=sqrtN_vals(i_N);
                    L=L_vals(i_N);
    
                    tlog=averaging_times{i_N,i_T};
                    t=tlog;
%                     tlin=averaging_times_lintime{i_N,i_T_lintime};
    
    
                    q_vals=qbin{i_N,i_T};
%                     q_vals_lintime=qbin_lintime{i_N,i_T_lintime};
                    i_q = find(q - q_vals < 1e-4,1);
                    
    
                    dispname=sprintf('$N = (%d)^2$', sqrtN);
                    q_indices = (i_q):length(q_vals):length(q_vals)*length(tlog);
%                     q_indices_lintime = (i_q):length(q_vals_lintime):length(q_vals_lintime)*length(tlin);
%                     gamma_cur = param_TCFSpin_omega_1_q_DO{i_N,i_T_lintime}(2*i_q+1);
%                     omega_1_cur = param_TCFSpin_omega_1_q_DO{i_N,i_T_lintime}(2*i_q+2);
%                     omega_0_cur = sqrt(omega_1_cur^2 + .25 * gamma_cur^2);
                    
                    
    
%                     t=averaging_times{i_N,i_T};

                    if (plot_type_cur == "m")
                        cf = gxx{i_N,i_T}(q_indices) + gyy{i_N,i_T}(q_indices);
%                         cf_lintime = gxx_lintime{i_N,i_T_lintime}(q_indices_lintime) ...
%                             + gyy_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "mpar")
                        cf = gmparmpar{i_N,i_T}(q_indices);
%                         cf_lintime = gmparmpar_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "mperp")
                        cf = gmperpmperp{i_N,i_T}(q_indices);
%                         cf_lintime = gmperpmperp_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "te")
                        cf = gtt{i_N,i_T}(q_indices);
%                         cf_lintime = gtt_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "w")
                        cf = gww{i_N,i_T}(q_indices);
%                         cf_lintime = gww_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    end
%                     [t,cf] = combine_cf(tlog,cf,tlin,cf_lintime);
                    coeffs_omega_1_cur=fit_DampedOscillator_RealSpace(t,real(cf),10,1,'omega_1');
                    gamma_cur = coeffs_omega_1_cur(1);
                    omega_1_cur = coeffs_omega_1_cur(2);
                    omega_0_cur = sqrt(omega_1_cur^2 + .25 * gamma_cur^2);
                    c_cur = omega_0_cur / q;
                    
%                     tau=res_factor/gamma_cur;
                    tau=res_factor/omega_1_cur;

%                     t=t(1:100);
%                     cf=cf(1:100);
                    resolution_vals=res_function(t,tau);
                    
%                     cf=fitfunc_DO(t,1,[1e-3,1e-2]);
                    [ft_vals,om_vals]=FT_correlation(t, cf .* resolution_vals, 1e6);
                    % Rough estimate of the integral
                    S_q = sum(real(ft_vals(1:end))) * (om_vals(2) - om_vals(1));
                    
%                     ft_vals = smoothen(abs(real(ft_vals)),1,3);
        %             S_q = chimxq_av{i_N,i_T}(i_q) + chimyq_av{i_N,i_T}(i_q);
    
                    loglog(om_vals / (c_cur * q ), c_cur * q * abs(real(ft_vals)) / S_q,...
                        'Color',c_map(i_N - i_N_start + 1,:),'DisplayName',dispname, ...
                        'LineStyle', '-', 'LineWidth',2, ...
                        'Marker', '.', 'MarkerSize', 8);
                    hold on;
                end
                xlim_vec=[3e-1,4];
                xlim(xlim_vec);
                
                ft_vals=fitfunc_DO_reciprocal(om_vals,1,[gamma_cur,omega_0_cur]);
                S_q = sum(real(ft_vals(1:end))) * (om_vals(2) - om_vals(1));
                ft_scaled=c_cur * q * real(ft_vals) / S_q;
                loglog(om_vals / (c_cur * q ), ft_scaled,...
                        'Color','black', 'DisplayName','DO Fit', ... 'Color',c_map(i_N - i_N_start + 2,:),
                        'LineStyle', '-.', 'LineWidth',2, ...
                        'Marker', 'none', 'MarkerSize', 5);
                
                

%                 om_vals_c_q=linspace(2,6);
%                 loglog(om_vals_c_q, om_vals_c_q.^(-(3-eta)),...
%                         'Color','Black','DisplayName','$\sim \omega^{-(3-\eta)}$', ...
%                         'LineStyle', '--', 'LineWidth',2, ...
%                         'Marker', 'none', 'MarkerSize', 5);

                y_vals=logspace(log10(xlim_vec(1)),log10(xlim_vec(2)),1e2);
%                 NF_Psi = c_cur * q * NelsonFisher_Psi(y_vals,eta,n_trunc) / S_q;
                NF_Psi = c_cur * q * NelsonFisher_Psi(y_vals,eta,n_trunc);
%                 NF_Psi = NelsonFisher_Psi(y_vals,eta,n_trunc);
%                 NF_Psi = max(ft_scaled)/max(NF_Psi)*NelsonFisher_Psi(y_vals,eta,n_trunc);
%                 loglog(om_vals,2*pi^2*eta^2*om_vals.^(-3-eta),...
                loglog(y_vals,NF_Psi,'--',...
                        'Color','Black','DisplayName','NF', ...
                        'LineWidth',2, ... 'Color',c_map(i_N - i_N_start + 3,:), 
                        'Marker', 'none', 'MarkerSize', 5);

                y_max=4*max([max(ft_scaled),NF_Psi]);
                ylim([1e-2/( max(om_vals) / (c_cur * q)), y_max]);

                h_axis = gca;
                hXLabel = xlabel('$\omega / cq$','interpreter','latex');
                if (plot_type_cur == "m")
                    hYLabel = ylabel('$ cq S_{m}(q,\omega) / S_{m}(q)$','interpreter','latex');
                elseif (plot_type_cur == "mpar")
                    hYLabel = ylabel('$ cq S_{m\parallel}(q,\omega) / S_{m\parallel}(q)$','interpreter','latex');
                elseif (plot_type_cur == "mperp")
                    hYLabel = ylabel('$ cq S_{m\perp}(q,\omega) / S_{m\perp}(q)$','interpreter','latex');
                elseif (plot_type_cur == "te")
                    hYLabel = ylabel('$ cq S_{\theta}(q,\omega) / S_{\theta}(q)$','interpreter','latex');
                elseif (plot_type_cur == "w")
                    hYLabel = ylabel('$ cq S_{w}(q,\omega) / S_{w}(q)$','interpreter','latex');
                end
                
    
    
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Legend, axes etc
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hLegend = legend('Location', 'northwest','interpreter','latex','FontSize', fontsize_annotation,...
                    'NumColumns',1);
    
                set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
                    'XMinorTick', 'on', 'YMinorTick', 'on', ...
                    'XGrid', 'off', 'YGrid', 'off', ...
                    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
                    'XTick', [0:.5:1,2:10], 'YTick', 10.^[-10:5], ...
                    'LineWidth', .5)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Font
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Annotations
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                annotation_str = {sprintf('$T = %.3f$',T),sprintf('$q = %.3f$',q),...
                    sprintf('$\\eta = %.2f$',eta)}; % ,sprintf('$z = %.2f$',z)
                dim=[.68 .68 .1 .2];
                annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                    'interpreter','latex',...
                    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
                    'Color','black','FontSize', fontsize_annotation,...
                    'BackgroundColor','white');

                

                set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    
                figname=sprintf('%s/HighFreq/%s/%s_S_FS_q_%.3f_z_%.2f_T_%.3f_Nstart_%d_resfac_%d',basedir,plot_type_cur,curmodel,q,z,T,i_N_start,res_factor);
    %                 figname=sprintf('%s/SbyL_vsomL/%s_S_%s_FS_n_q_%d_z_%.2f_tau_%de3_T_%.3f',basedir,curmodel,plot_type_cur,n_q,z,tau_laplace/1e3,T);
                fprintf('Creating figure %s\n',figname)
                if(saveswitch == 1)
                   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
                end
            end
        end
        

        %% Section 4 Plot to study the Nelson-Fisher high-frequency with c, scaled to omega^{3-eta}
        c_map = linspecer(numel(sqrtN_vals));
        if (ismember(4,section_select))
            i_N_start = 3;
            q = qbin{i_N_start,1}(1);
            for i_T = T_select
                figure;
                T=T_vals(i_T);
                i_T_lintime=find_T_index(lintime_dataset_name,T);
                eta=eta_vals(i_T_lintime);
                c_cur = c_vals(i_T_lintime);
                c_map = linspecer(numel(i_N_start : 5));
                for i_N = i_N_start:numel(sqrtN_vals)
                    sqrtN=sqrtN_vals(i_N);
                    L=L_vals(i_N);
    
                    tlog=averaging_times{i_N,i_T};
                    tlin=averaging_times_lintime{i_N,i_T_lintime};
    
    
                    q_vals=qbin{i_N,i_T};
                    q_vals_lintime=qbin_lintime{i_N,i_T_lintime};
                    i_q = find(q - q_vals < 1e-4,1);
                    
    
                    dispname=sprintf('$N = (%d)^2$', sqrtN);
                    q_indices = (i_q):length(q_vals):length(q_vals)*length(tlog);
                    q_indices_lintime = (i_q):length(q_vals_lintime):length(q_vals_lintime)*length(tlin);
%                     gamma_cur = param_TCFSpin_omega_1_q_DO{i_N,i_T_lintime}(2*i_q+1);
%                     omega_1_cur = param_TCFSpin_omega_1_q_DO{i_N,i_T_lintime}(2*i_q+2);
%                     omega_0_cur = sqrt(omega_1_cur^2 + .25 * gamma_cur^2);
                    
                    tau=res_factor/gamma_cur;
    
%                     t=averaging_times{i_N,i_T};

                    if (plot_type_cur == "m")
                        cf = gxx{i_N,i_T}(q_indices) + gyy{i_N,i_T}(q_indices);
                        cf_lintime = gxx_lintime{i_N,i_T_lintime}(q_indices_lintime) ...
                            + gyy_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "mpar")
                        cf = gmparmpar{i_N,i_T}(q_indices);
                        cf_lintime = gmparmpar_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "mperp")
                        cf = gmperpmperp{i_N,i_T}(q_indices);
                        cf_lintime = gmperpmperp_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "te")
                        cf = gtt{i_N,i_T}(q_indices);
                        cf_lintime = gtt_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "w")
                        cf = gww{i_N,i_T}(q_indices);
                        cf_lintime = gww_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    end
                    [t,cf] = combine_cf(tlog,cf,tlin,cf_lintime);
                    coeffs_omega_1_cur=fit_DampedOscillator_RealSpace(t,real(cf),10,1,'omega_1');
                    gamma_cur = coeffs_omega_1_cur(1);
                    omega_1_cur = coeffs_omega_1_cur(1);
                    omega_0_cur = sqrt(omega_1_cur^2 + .25 * gamma_cur^2);
%                     coeffs_omega_1_cur(2*(i_q-1)+(1:2))=fit_DampedOscillator_RealSpace(t,real(cf),10,1,'omega_0');
%                     t=t(1:100);
%                     cf=cf(1:100);
                    resolution_vals=res_function(t,tau);
                    
%                     cf=fitfunc_DO(t,1,[1e-3,1e-2]);
                    [ft_vals,om_vals]=FT_correlation(t, cf .* resolution_vals, 1e6+1);
                    % Rough estimate of the integral
                    S_q = sum(real(ft_vals(1:end))) * (om_vals(2) - om_vals(1));
        %             S_q = chimxq_av{i_N,i_T}(i_q) + chimyq_av{i_N,i_T}(i_q);
    
                    loglog(om_vals / (c_cur * q ), ...
                        c_cur * q * abs(real(ft_vals)) / S_q .* (om_vals / (c_cur * q )).^(3-eta),...
                        'Color',c_map(i_N - i_N_start + 1,:),'DisplayName',dispname, ...
                        'LineStyle', 'none', 'LineWidth',2, ...
                        'Marker', '.', 'MarkerSize', 5);
                    hold on;
                end
                xlim([3e-1,1e2]);
                ylim([1e-2/( max(om_vals) / (c_cur * q)), Inf])

                h_axis = gca;
                hXLabel = xlabel('$\omega / cq$','interpreter','latex');
                if (plot_type_cur == "m")
                    hYLabel = ylabel('$ cq S_{m}(q,\omega) / S_{m}(q) (\omega/cq)^{3-\eta}$','interpreter','latex');
                elseif (plot_type_cur == "mpar")
                    hYLabel = ylabel('$ cq S_{m\parallel}(q,\omega) / S_{m\parallel}(q) (\omega/cq)^{3-\eta}$','interpreter','latex');
                elseif (plot_type_cur == "mperp")
                    hYLabel = ylabel('$ cq S_{m\perp}(q,\omega) / S_{m\perp}(q) (\omega/cq)^{3-\eta}$','interpreter','latex');
                elseif (plot_type_cur == "te")
                    hYLabel = ylabel('$ cq S_{\theta}(q,\omega) / S_{\theta}(q) (\omega/cq)^{3-\eta}$','interpreter','latex');
                elseif (plot_type_cur == "w")
                    hYLabel = ylabel('$ cq S_{w}(q,\omega) / S_{w}(q) (\omega/cq)^{3-\eta}$','interpreter','latex');
                end
                
    
    
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Legend, axes etc
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hLegend = legend('Location', 'northwest','interpreter','latex','FontSize', fontsize_annotation,...
                    'NumColumns',1);
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Font
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)
    
                annotation_str = {sprintf('$T = %.3f$',T),sprintf('$q = %.3f$',q),...
                    sprintf('$\\eta = %.2f$',eta),sprintf('$z = %.2f$',z)};
                dim=[.44 .68 .1 .2];
                annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                    'interpreter','latex',...
                    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
                    'Color','black','FontSize', fontsize_annotation,...
                    'BackgroundColor','white');
    
                figname=sprintf('%s/HighFreq/%s/%s_S_FS_omscaled_q_%.3f_z_%.2f_T_%.3f_Nstart_%d_resfac_%d',basedir,plot_type_cur,curmodel,q,z,T,i_N_start,res_factor);

                set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    %                 figname=sprintf('%s/SbyL_vsomL/%s_S_%s_FS_n_q_%d_z_%.2f_tau_%de3_T_%.3f',basedir,curmodel,plot_type_cur,n_q,z,tau_laplace/1e3,T);
                fprintf('Creating figure %s\n',figname)
                if(saveswitch == 1)
                   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
                end
            end
        end
       
        %% Section 5 Same as Section 4, but with 2eta instead of eta
        c_map = linspecer(numel(sqrtN_vals));
        if (ismember(5,section_select))
            i_N_start = 3;
            q = qbin{i_N_start,1}(1);
            for i_T = T_select
                figure;
                T=T_vals(i_T);
                i_T_lintime=find_T_index(lintime_dataset_name,T);
                eta=eta_vals(i_T_lintime);
                c_cur = c_vals(i_T_lintime);
                c_map = linspecer(numel(i_N_start : 5));
                for i_N = i_N_start:numel(sqrtN_vals)
                    sqrtN=sqrtN_vals(i_N);
                    L=L_vals(i_N);
    
                    tlog=averaging_times{i_N,i_T};
                    tlin=averaging_times_lintime{i_N,i_T_lintime};
    
    
                    q_vals=qbin{i_N,i_T};
                    q_vals_lintime=qbin_lintime{i_N,i_T_lintime};
                    i_q = find(q - q_vals < 1e-4,1);
                    
    
                    dispname=sprintf('$N = (%d)^2$', sqrtN);
                    q_indices = (i_q):length(q_vals):length(q_vals)*length(tlog);
                    q_indices_lintime = (i_q):length(q_vals_lintime):length(q_vals_lintime)*length(tlin);
                    gamma_cur = param_TCFSpin_omega_1_q_DO{i_N,i_T_lintime}(2*i_q+1);
                    omega_1_cur = param_TCFSpin_omega_1_q_DO{i_N,i_T_lintime}(2*i_q+2);
                    omega_0_cur = sqrt(omega_1_cur^2 + .25 * gamma_cur^2);
                    
                    tau=res_factor/gamma_cur;
    
%                     t=averaging_times{i_N,i_T};

                    if (plot_type_cur == "m")
                        cf = gxx{i_N,i_T}(q_indices) + gyy{i_N,i_T}(q_indices);
                        cf_lintime = gxx_lintime{i_N,i_T_lintime}(q_indices_lintime) ...
                            + gyy_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "mpar")
                        cf = gmparmpar{i_N,i_T}(q_indices);
                        cf_lintime = gmparmpar_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "mperp")
                        cf = gmperpmperp{i_N,i_T}(q_indices);
                        cf_lintime = gmperpmperp_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "te")
                        cf = gtt{i_N,i_T}(q_indices);
                        cf_lintime = gtt_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    elseif (plot_type_cur == "w")
                        cf = gww{i_N,i_T}(q_indices);
                        cf_lintime = gww_lintime{i_N,i_T_lintime}(q_indices_lintime);
                    end
                    [t,cf] = combine_cf(tlog,cf,tlin,cf_lintime);

%                     t=t(1:100);
%                     cf=cf(1:100);
                    resolution_vals=res_function(t,tau);
                    
%                     cf=fitfunc_DO(t,1,[1e-3,1e-2]);
                    [ft_vals,om_vals]=FT_correlation(t, cf .* resolution_vals, 1e6+1);
                    % Rough estimate of the integral
                    S_q = sum(real(ft_vals(1:end))) * (om_vals(2) - om_vals(1));
        %             S_q = chimxq_av{i_N,i_T}(i_q) + chimyq_av{i_N,i_T}(i_q);
    
                    loglog(om_vals / (c_cur * q ), ...
                        c_cur * q * abs(real(ft_vals)) / S_q .* (om_vals / (c_cur * q )).^(3-2*eta),...
                        'Color',c_map(i_N - i_N_start + 1,:),'DisplayName',dispname, ...
                        'LineStyle', 'none', 'LineWidth',2, ...
                        'Marker', '.', 'MarkerSize', 5);
                    hold on;
                end
                xlim([3e-1,1e2]);
                ylim([1e-2/( max(om_vals) / (c_cur * q)), Inf])

                h_axis = gca;
                hXLabel = xlabel('$\omega / cq$','interpreter','latex');
                if (plot_type_cur == "m")
                    hYLabel = ylabel('$ cq S_{m}(q,\omega) / S_{m}(q) (\omega/cq)^{3-2\eta}$','interpreter','latex');
                elseif (plot_type_cur == "mpar")
                    hYLabel = ylabel('$ cq S_{m\parallel}(q,\omega) / S_{m\parallel}(q) (\omega/cq)^{3-2\eta}$','interpreter','latex');
                elseif (plot_type_cur == "mperp")
                    hYLabel = ylabel('$ cq S_{m\perp}(q,\omega) / S_{m\perp}(q) (\omega/cq)^{3-2\eta}$','interpreter','latex');
                elseif (plot_type_cur == "te")
                    hYLabel = ylabel('$ cq S_{\theta}(q,\omega) / S_{\theta}(q) (\omega/cq)^{3-2\eta}$','interpreter','latex');
                elseif (plot_type_cur == "w")
                    hYLabel = ylabel('$ cq S_{w}(q,\omega) / S_{w}(q) (\omega/cq)^{3-2\eta}$','interpreter','latex');
                end
                
    
    
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Legend, axes etc
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hLegend = legend('Location', 'northwest','interpreter','latex','FontSize', fontsize_annotation,...
                    'NumColumns',1);
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Font
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)
    
                annotation_str = {sprintf('$T = %.3f$',T),sprintf('$q = %.3f$',q),...
                    sprintf('$\\eta = %.2f$',eta),sprintf('$z = %.2f$',z)};
                dim=[.44 .68 .1 .2];
                annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                    'interpreter','latex',...
                    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
                    'Color','black','FontSize', fontsize_annotation,...
                    'BackgroundColor','white');

                set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    
                figname=sprintf('%s/HighFreq/%s/%s_S_FS_omscaled_2eta_q_%.3f_z_%.2f_T_%.3f_Nstart_%d_resfac_%d',basedir,plot_type_cur,curmodel,q,z,T,i_N_start,res_factor);
    %                 figname=sprintf('%s/SbyL_vsomL/%s_S_%s_FS_n_q_%d_z_%.2f_tau_%de3_T_%.3f',basedir,curmodel,plot_type_cur,n_q,z,tau_laplace/1e3,T);
                fprintf('Creating figure %s\n',figname)
                if(saveswitch == 1)
                   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
                end
            end
        end
        
    end
end
   