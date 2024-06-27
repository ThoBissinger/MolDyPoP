clear
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/TCF/FFT_FineStructure',fig_base);

section_select = 1;

for i_model = [1]
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";

        data=matfile(mxydata_AdjustedTime_name);
        
        fitdata=matfile(mxyfit_name);
        fitdata_FS=matfile(mxyfit_FS_name);
        fitdata_TCF=matfile(mxyfit_TCF_q_name);

        lintime_dataset_name='dynamics_mxy_LinearTime';

        i_N=5;
        i_T=4;
        q_select=[1,3,9,15];
        xlim_vec=[0 3];
        
        L_vals=[9.25,18.5,37,74,148];

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
        pairs=[1,3,8; 2,3,8; 3,3,8; 4,3,8]; % T=.91
        xlim_vec=[0 2];
    end
    res_factor=3;
    res_function = resolution_Gauss;
    res_function=@(t,tau) resolution_Laplace_pleateau(t,tau,tau);

    i_N_start = 2;
    i_N_max = 4;
    n_trunc = 4; % For NF Psi
%     plot_type=["w", "mperp", "te","mpar", "m"];
%     plot_type=["mperp", "te","mpar", "m"];
    plot_type=["mperp", "m"];
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
    for i_plot_type = 1:numel(plot_type)
        plot_type_cur = plot_type(i_plot_type);
        %% Section 1 Plot to study the Nelson-Fisher high-frequency behavior
        N_q=numel(q_select);
        c_map = linspecer(N_q);
        figure
        hold on;
        max_ft_vec=zeros(1,N_q);
        S_q_vec=zeros(1,N_q);
        om_max_vec=zeros(1,N_q);
        q_vals = qbin{i_N,i_T};
        T = T_vals(i_T);
        sqrtN = sqrtN_vals(i_N);
        L = L_vals(i_N);
        t=averaging_times{i_N,i_T};
        N_t = numel(t);
        N_q = numel(q_vals);
        
        for ind_q = 1:N_q
            i_q = q_select(ind_q);
            q = q_vals(i_q);
            
            q_indices=i_q:N_q:N_q*N_t;    
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
            resolution_vals=res_function(t,tau);
                
            [ft_vals,om_vals]=FT_correlation(t, cf .* resolution_vals, 1e5);
            ft_vals=real(ft_vals);
            S_q = cf(1);
            [ft_max,i_ft_max] = max(ft_vals);
            
            max_ft_vec(i_pair) = ft_max;
            S_q_vec(i_pair) = S_q;
            om_max_vec(i_pair) = om_vals(i_ft_max);
    
            dispname = sprintf('$N = (%d)^2$',sqrtN);
            loglog(om_vals, real(ft_vals)/ L^z / S_q,...
                'Color',c_map(i_pair,:),'DisplayName',dispname, ...
                'LineStyle', '-', 'LineWidth',2, ...
                'Marker', 'none', 'MarkerSize', 3);
            hold on;
        end
        absM_vec=cell2mat(absM_av(:,i_T));
        eta_fitob=fit_eta_Magnetization_FS(absM_vec,L_vals);
        eta=eta_fitob.eta;
%         y_vals=logspace(log10(xlim_vec(1)),log10(xlim_vec(2)),1e2);
%         NF_Psi = NelsonFisher_Psi(y_vals,eta,n_trunc);
% %                 loglog(om_vals,2*pi^2*eta^2*om_vals.^(-3-eta),...
%         loglog(y_vals,NF_Psi,...
%                 'Color','Black','DisplayName','NF', ...
%                 'LineStyle', '--', 'LineWidth',2, ...
%                 'Marker', 'none', 'MarkerSize', 5);

        xlim(xlim_vec);
        ylim([0 Inf]);
%         ylim([.3*min(NF_Psi),10*max(NF_Psi)]);
        h_axis = gca;
        hXLabel = xlabel('$\omega $','interpreter','latex');
        if (plot_type_cur == "m")
            hYLabel = ylabel('$S_{m}(q,\omega) / \chi_{m}(q) L^z$','interpreter','latex');
        elseif (plot_type_cur == "mpar")
            hYLabel = ylabel('$S_{m\parallel}(q,\omega) / \chi_{m\parallel}(q) L^z$','interpreter','latex');
        elseif (plot_type_cur == "mperp")
            hYLabel = ylabel('$S_{m\perp}(q,\omega) / \chi_{m\perp}(q) L^z$','interpreter','latex');
        elseif (plot_type_cur == "te")
            hYLabel = ylabel('$S_{\theta}(q,\omega) / \chi_{\theta}(q) L^z$','interpreter','latex');
        elseif (plot_type_cur == "w")
            hYLabel = ylabel('$S_{w}(q,\omega) / \chi_{w}(q) L^z$','interpreter','latex');
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
        dim=[.72 .2 .1 .2];
        annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
            'interpreter','latex',...
            'VerticalAlignment','middle', 'HorizontalAlignment','left',...
            'Color','black','FontSize', fontsize_annotation,...
            'BackgroundColor','white');

        set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

%                 figname=sprintf('%s/HighFreq/%s/%s_S_FS_n_q_%d_z_%.2f_T_%.3f_Nstart_%d_resfac_%d',basedir,plot_type_cur,curmodel,n_q,z,T,i_N_start,res_factor);
        figname=sprintf('%s/%s_%s_TCF_FFT_FS_T_%.3f_nq_%d',basedir,plot_type_cur,curmodel,T,N_q);
%                 figname=sprintf('%s/SbyL_vsomL/%s_S_%s_FS_n_q_%d_z_%.2f_tau_%de3_T_%.3f',basedir,curmodel,plot_type_cur,n_q,z,tau_laplace/1e3,T);
        fprintf('Creating figure %s\n',figname)
        if(saveswitch == 1)
           exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');

           set(gca,'YScale','log')
           figname=sprintf('%s/%s_%s_TCF_FFT_FS_logy_T_%.3f_nq_%d',basedir,plot_type_cur,curmodel,T,N_q);
           ylim([1e-2 Inf])
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end

        
    end
end
