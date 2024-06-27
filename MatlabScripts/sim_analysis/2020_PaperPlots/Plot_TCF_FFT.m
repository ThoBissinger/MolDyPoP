%% 0 Initialization

run initialization_script;
saveswitch=0;
fontsize_axis=15;
fontsize_annotation=15;
basedir=sprintf('%s/plots/FFT',fig_base);
fig_select=[1:3];
for i_model = [1]

    % figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=matfile(mxydata_name);
        fitdata=matfile(mxyfit_TCF_q_name);
        fit_FS=matfile(mxyfit_FS_name);
        
        T_select=[2,4,8,10,11:17];
        xmax = .45;
        L_vals = [9.25,18.5,37,74,148];
        z=1;
        n_q = 1;
        
        res_factor=10;
        tau_laplace=4e3;
    elseif (i_model == 2)
        curmodel="xy";
        curtitle="SXY model";
        
        data=matfile(xydata_name);
        fitdata=matfile(xyfit_TCF_q_name);
        fit_FS=matfile(xyfit_FS_name);
        
        T_select=[2,4,8,10,11:19];
        xmax=2.5 ;
        L_vals=[16,32,64,128,256];
        z=1;
        n_q = 1;
        
        res_factor=10;
        tau_laplace=4e3;
        
        
        
        
    elseif (i_model == 3)
        curmodel="fmxy";
        curtitle="FMXY model";
        
        data=matfile(fmxydata_name);
        fitdata=matfile(fmxyfit_TCF_q_name);
        fit_FS=matfile(fmxyfit_FS_name);
        
        T_select=[2,4,8,10,11:17];
        xmax = .45;
        L_vals = [9.25,18.5,37,74,148];
        z=1;
        n_q = 1;
        
        res_factor=10;
        tau_laplace=4e3;
    end

    plot_type=["w", "mperp", "te","mpar", "m"];
    plot_type=["mperp", "te","mpar", "m"];
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

    eta_vals=fit_FS.('eta_vals');
    sqrtN_vals=data.('sqrtN_vals');
    T_vals=data.('T_vals');

    param_TCFSpin_omega_1_q_DO = fitdata.('param_TCFSpin_omega_1_q_DO');

    for i_plot_type = 1:numel(plot_type)
        plot_type_cur = plot_type(i_plot_type);
        
        %% 1 S(q,omega) q^{3-eta} vs omega/q
        if(ismember(1,fig_select))
    
            for i_T = T_select
%                 figure(i_T + 1e4*i_model + 1e5 * i_plot_type);
                figure;
                i_N=5;
    
                sqrtN=sqrtN_vals(i_N);
                T=T_vals(i_T);
                eta=eta_vals(i_T);
                
                t=averaging_times{i_N,i_T};
    
                q_vals=qbin{i_N,i_T};
                q_select = 1:6;
                c_map = linspecer(numel(q_select));
                for index_q = 1:numel(q_select)
                    i_q = q_select(index_q);
                    q = q_vals(i_q );
                    dispname=sprintf('$q = %.3f$', q);
                    q_indices = (i_q):length(q_vals):length(q_vals)*length(t);
                    
                    gamma_cur=param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+1);
                    omega_1_cur=param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+2);
                    
                    tau=res_factor/gamma_cur;

                    
                    resolution_vals=resolution_Laplace(t,tau);
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
                    [ft_vals,om_vals]=FT_correlation(t, cf .* resolution_vals, 0);
                    plot(om_vals/q,abs(ft_vals)*q^(3-eta),...
                        'Color',c_map(index_q,:),'DisplayName',dispname, ...
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
                hLegend = legend('Location', 'Northeast','interpreter','latex','FontSize', fontsize_annotation,...
                    'NumColumns',1);
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Font
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)
    
    
                annotation_str = {sprintf('$N = (%d)^2$',sqrtN),sprintf('$T = %.3f$',T),sprintf('$\\eta = %.2f$',eta)};
                dim=[.16 .7 .1 .2];
                annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                    'interpreter','latex',...
                    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
                    'Color','Red','FontSize', fontsize_annotation,...
                    'BackgroundColor','white');
    
                figname=sprintf('%s/S_by_om_q_scaling/%s/%s_S_%s_sqrtN_%d_T_%.3f_resfac_%d',basedir,plot_type_cur,curmodel,plot_type_cur,sqrtN,T,res_factor);
                fprintf('Creating figure %s\n',figname)
                if(saveswitch == 1)
                   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
                end
            end
        end
    
    
    
        %% 2 Figure S(q,om)/S(q) vs om
        if(ismember(2,fig_select))
    
            i_N_start = 3;
            q = qbin{i_N_start,1}(1);
            for i_T = T_select
%                 figure(i_T + 100 + 1e4*i_model + 1e5 * i_plot_type);
                figure;
                T=T_vals(i_T);
                eta=eta_vals(i_T);
                
                c_map = linspecer(numel(i_N_start : 5));
                for i_N = i_N_start:5
                    sqrtN=sqrtN_vals(i_N);
                    L=L_vals(i_N);
                    
                    t=averaging_times{i_N,i_T};
    
                    q_vals=qbin{i_N,i_T};
                    i_q = find(min(abs(q_vals - q)) == abs(q_vals - q));
    
                    dispname=sprintf('$N = (%d)^2$', sqrtN);
                    q_indices = (i_q):length(q_vals):length(q_vals)*length(t);
                    
                    gamma_cur=param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+1);
                    omega_1_cur=param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+2);
                    
                    tau=res_factor/gamma_cur;

                    t=averaging_times{i_N,i_T};
                    resolution_vals=resolution_Laplace(t,tau);
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
                    [ft_vals,om_vals]=FT_correlation(t, cf .* resolution_vals, 0);
                    % Rough estimate of the integral
                    S_q = sum(abs(ft_vals(1:end))) * (om_vals(2) - om_vals(1));
        %             S_q = chimxq_av{i_N,i_T}(i_q) + chimyq_av{i_N,i_T}(i_q);
    
                    plot(om_vals,abs(ft_vals) / S_q,...
                        'Color',c_map(i_N - i_N_start + 1,:),'DisplayName',dispname, ...
                        'LineStyle', '-', 'LineWidth',2, ...
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
                hLegend = legend('Location', 'Northeast','interpreter','latex','FontSize', fontsize_annotation,...
                    'NumColumns',1);
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Font
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)
    
                annotation_str = {sprintf('$T = %.3f$',T),sprintf('$q = %.3f$',q),...
                    sprintf('$\\eta = %.2f$',eta)};
                dim=[.68 .45 .1 .2];
                annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                    'interpreter','latex',...
                    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
                    'Color','Black','FontSize', fontsize_annotation,...
                    'BackgroundColor','white');
    
                figname=sprintf('%s/S_vs_om_fixed_q/%s/%s_S_%s_FS_T_%.3f_Nstart_%d_resfac_%d',basedir,plot_type_cur,curmodel,plot_type_cur,T,i_N_start,res_factor);
%                 figname=sprintf('%s/S_vs_om_fixed_q/%s_S_%s_FS_T_%.3f',basedir,curmodel,plot_type_cur,T);
                fprintf('Creating figure %s\n',figname)
                if(saveswitch == 1)
                   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
                end
            end
        end
    
    
    %% 3 Figure S(q,om)/L^zS(q) vs om L^z
        if(ismember(3,fig_select))
    
            i_N_start = 3;
            q = qbin{i_N_start,1}(n_q);
            for i_T = T_select
%                 figure(i_T + 300 + 1e4*i_model + 1e5 * i_plot_type);
                figure;
                T=T_vals(i_T);
                eta=eta_vals(i_T);
                c_map = linspecer(numel(i_N_start : 5));
                for i_N = i_N_start:5
                    sqrtN=sqrtN_vals(i_N);
                    L=L_vals(i_N);

                    t=averaging_times{i_N,i_T};
    

                    q_vals=qbin{i_N,i_T};
                    i_q = n_q;
    
                    dispname=sprintf('$N = (%d)^2$', sqrtN);
                    q_indices = (i_q):length(q_vals):length(q_vals)*length(t);
                    gamma_cur=param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+1);
                    omega_1_cur=param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+2);
                    
                    tau=res_factor/gamma_cur;

                    t=averaging_times{i_N,i_T};
                    resolution_vals=resolution_Laplace(t,tau);
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
                    [ft_vals,om_vals]=FT_correlation(t, cf .* resolution_vals, 1e6);
                    % Rough estimate of the integral
                    S_q = sum(abs(ft_vals(1:end))) * (om_vals(2) - om_vals(1));
        %             S_q = chimxq_av{i_N,i_T}(i_q) + chimyq_av{i_N,i_T}(i_q);
    
                    loglog(om_vals*L^z,abs(ft_vals)/ L^z / S_q,...
                        'Color',c_map(i_N - i_N_start + 1,:),'DisplayName',dispname, ...
                        'LineStyle', '-', 'LineWidth',2, ...
                        'Marker', 'none', 'MarkerSize', 5);
                    hold on;
                end
                xlim([.5,Inf]);
                loglog(logspace(1,3.5),.1*logspace(1,3.5).^(-1),...
                        'Color','Black','DisplayName','$(\omega L^z)^{-1}$', ...
                        'LineStyle', '--', 'LineWidth',2, ...
                        'Marker', 'none', 'MarkerSize', 5);
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
                hLegend = legend('Location', 'Northeast','interpreter','latex','FontSize', fontsize_annotation,...
                    'NumColumns',1);
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Font
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)
    
                annotation_str = {sprintf('$T = %.3f$',T),sprintf('$n_q = %d$',i_q),...
                    sprintf('$\\eta = %.2f$',eta),sprintf('$z = %.2f$',z)};
                dim=[.685 .4 .1 .2];
                annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                    'interpreter','latex',...
                    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
                    'Color','Red','FontSize', fontsize_annotation,...
                    'BackgroundColor','white');
    
                figname=sprintf('%s/SbyL_vsomL/%s/%s_S_%s_FS_n_q_%d_z_%.2f_T_%.3f_Nstart_%d_resfac_%d',basedir,plot_type_cur,curmodel,plot_type_cur,n_q,z,T,i_N_start,res_factor);
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