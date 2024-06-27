run initialization_script;
saveswitch=1;
fontsize_axis=15;
fontsize_annotation=15;
basedir=sprintf('%s/plots/TCF/mperp_TCF_FFT_Examples',fig_base);

for i_model = [1]

    % figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=matfile(mxydata_AdjustedTime_name);
        fitdata=matfile(mxyfit_TCF_q_name);

        % pairs of i_N, i_T, i_q
        qfixed_pairs=[4,2,1; 4,4,1; 4,6,1];
        Tfixed_pairs=[4,4,1; 4,4,2; 4,4,3];
        L_vals=[9.25,18.5,37,74,148];
        runmax=500;
    elseif (i_model == 2)
        curmodel="xy";
        curtitle="SXY model";
        
        data=matfile(xydata_name);
        fitdata=matfile(xyfit_TCF_q_name);
        
        T_select=[17:3:38];
        
        L_vals=[16,32,64,128,256];
        
        q_select=1:6;
        runmax=250;
    elseif (i_model == 3)
        curmodel="fmxy";
        curtitle="FMXY model";
        
        data=matfile(fmxydata_name);
        fitdata=matfile(fmxyfit_TCF_q_name);
        
        qfixed_pairs=[4,2,1; 4,4,1; 4,6,1];
        Tfixed_pairs=[4,4,1; 4,4,2; 4,4,3];
        
        L_vals=[9.25,18.5,37,74,148];
        runmax=500;
    elseif (i_model == 4)
        curmodel="xy_s";
        curtitle="XY S model";
        model_str="XY";

        data=matfile(xysdata_AdjustedTime_name);

        fitdata=matfile(xysfit_FSMag_name);
        fitdata_FS=matfile(xysfit_FSMag_name);
        fitdata_TCF=matfile(xysfit_TCF_q_name);
        
        lintime_dataset_name='dynamics_xy_s_AdjustedTime';

        
        qfixed_pairs=[4,1,1; 4,3,1; 4,6,1];
        Tfixed_pairs=[4,13,1; 4,13,2; 4,13,3];
        
        L_vals=[16, 32, 64, 128];
        runmax=500;
    end
    

    res_factor=4;
    res_function=@(t,tau) resolution_Laplace_pleateau(t,tau,tau);
    
    T_vals=data.('T_vals');
    sqrtN_vals=data.('sqrtN_vals');
    gmperpmperp=data.('gmperpmperp');
    TCF_times=data.('TCF_times');
    qbin = data.('qbin');
    
    param_TCFSpin_omega_1_q_DO = fitdata.('param_TCFSpin_omega_1_q_DO');
    
    i_N = numel(sqrtN_vals);
    q_vals = qbin{i_N,1};

%     [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% 1 FFT for various T at fixed q
    n_pairs=numel(qfixed_pairs(:,1));
    c_map = linspecer(n_pairs);
    h_cf_plot=cell(1,n_pairs);
    figure
    for i_pair = 1:n_pairs
        
        i_N = qfixed_pairs(i_pair,1);
        i_T = qfixed_pairs(i_pair,2);
        i_q = qfixed_pairs(i_pair,3);

        sqrtN = sqrtN_vals(i_N);
        T = T_vals(i_T);
        q_vals = qbin{i_N,i_T};
        q=q_vals(i_q);

        gamma_cur=param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+1);
        omega_cur=abs(param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+2));
        sigma=1/(res_factor*gamma_cur);
        
        t_vals=TCF_times{i_N,i_T};
        N_q=length(qbin{i_N,i_T});
        N_t=length(t_vals);
        q_indices = (i_q):N_q:N_t*N_q; 
        cf_vals=real(gmperpmperp{i_N,i_T}(q_indices));
    
        dispname=sprintf('$T = %.3f$',T);
        tau=res_factor/gamma_cur;
        resolution_vals=res_function(t,tau);
            
        [ft_vals,om_vals]=FT_correlation(t, cf .* resolution_vals, 1e5);
        h_cf_plot{i_pair}=plot(om_vals,ft_vals); 
        hold on; 
        set(h_cf_plot{i_pair},'LineStyle','-','LineWidth',2,...
            'DisplayName',dispname,'Color',c_map(i_pair,:));
    
        
    end
    xlim([0,2*omega_cur]);
        
    hXLabel = xlabel('$\omega$','interpreter','latex');
    hYLabel = ylabel('$S_{m\perp}(q,\omega)$','interpreter','latex');

    hlegend = legend('location','northeast','interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);
    % Font
    set([gca,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize',fontsize_axis)


    pos = get(gca, 'position');
    dim = [pos(3)-.11 .55*pos(2) pos(3) .75*pos(4)]; %top
    str = {sprintf('$N = (%d)^2$',sqrtN), sprintf('$q = %.3f$',q)}; %sprintf('$T = %.3f$',T), 
    h_ann=annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);
    
    if(saveswitch == 1)
        figname=sprintf('%s/%s_qfixed_sqrtN_%d_q_%.3f_resfactor_%.3f',basedir,curmodel,sqrtN,q,res_factor);
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    

    %% 2 FFT for various q at fixed T
    n_pairs=numel(Tfixed_pairs(:,1));
    c_map = linspecer(n_pairs);
    h_cf_plot=cell(size(Tfixed_pairs));
    figure(200*i_model + 1)
    for i_pair = 1:n_pairs
        
        i_N = Tfixed_pairs(i_pair,1);
        i_T = Tfixed_pairs(i_pair,2);
        i_q = Tfixed_pairs(i_pair,3);

        sqrtN = sqrtN_vals(i_N);
        T = T_vals(i_T);
        q_vals = qbin{i_N,i_T};
        q=q_vals(i_q);

        gamma_cur=param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+1);
        omega_cur=abs(param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+2));
        sigma=1/(res_factor*gamma_cur);
        
        t_vals=TCF_times{i_N,i_T};
        N_q=length(qbin{i_N,i_T});
        N_t=length(t_vals);
        q_indices = (i_q):N_q:N_t*N_q; 
        cf_vals=real(gmperpmperp{i_N,i_T}(q_indices));
    
        dispname=sprintf('$q = %.3f$',q);
        tau=res_factor/gamma_cur;
        resolution_vals=res_function(t,tau);
            
        [ft_vals,om_vals]=FT_correlation(t, cf .* resolution_vals, 1e5);
%         [ft_vals,om_vals]=FT_correlation(t_vals,cf_vals.*exp(-.5*t_vals.^2/sigma^2), 1e6);
        h_cf_plot{i_pair}=plot(om_vals,ft_vals); hold on; 
        set(h_cf_plot{i_pair},'LineStyle','-','LineWidth',2,...
            'DisplayName',dispname,'Color',c_map(i_pair,:));
        hold on;
    
        
    end
    xlim([0,1.5*omega_cur]);
        
    hXLabel = xlabel('$\omega$','interpreter','latex');
    hYLabel = ylabel('$S_{m\perp}(q,\omega)$','interpreter','latex');

    hlegend = legend('location','northeast','interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);
    % Font
    set([gca,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize',fontsize_axis)


    pos = get(gca, 'position');
    dim = [pos(3)-.11 .55*pos(2) pos(3) .75*pos(4)]; %top
    str = {sprintf('$N = (%d)^2$',sqrtN), sprintf('$T = %.3f$',T)}; %sprintf('$T = %.3f$',T), 
    h_ann=annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);
    
    if(saveswitch == 1)
        figname=sprintf('%s/%s_Tfixed_sqrtN_%d_T_%.3f_resfactor_%.3f',basedir,curmodel,sqrtN,T,res_factor);
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    


end
    