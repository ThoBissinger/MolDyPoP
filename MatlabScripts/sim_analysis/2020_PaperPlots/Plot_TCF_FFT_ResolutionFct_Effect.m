run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/TCF/mperp_TCF_Resolution',fig_base);

for i_model = [1]

    % figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=matfile(mxydata_AdjustedTime_name);
        fitdata=matfile(mxyfit_TCF_q_name);
        
        % pairs of i_N, i_T, i_q
        pairs=[4,3,1; 4,3,4; 4,3,7];
%         pairs=[4,10,1; 3,5,1; 5,10,3;...
%             2,10,3; 2,10,1; 3,10,1];
        res_vals=[.01,.05,.1,.3,.6,1,2];
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
        
        data=matfile(fmxydata_AdjustedTime_name);
        fitdata=matfile(fmxyfit_TCF_q_name);
        
        pairs=[4,10,3; 3,10,1; 5,10,6];

        L_vals=[9.25,18.5,37,74,148];
        runmax=500;
    end
    

    
    T_vals=data.('T_vals');
    sqrtN_vals=data.('sqrtN_vals');
    gmperpmperp=data.('gmperpmperp');
    TCF_times=data.('TCF_times');
    qbin = data.('qbin');
    
    param_TCFSpin_omega_1_q_DO = fitdata.('param_TCFSpin_omega_1_q_DO');
    
    fitfunc_DO=@(times,corrfunc_1,c) corrfunc_1 * exp(-c(1) * times/2) .* (cos(c(2) * times) + .5*c(1)/c(2)*sin(c(2) * times));  
    
    i_N = numel(sqrtN_vals);
    q_vals = qbin{i_N,1};

%     [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% 1 Behavior in time domain, Laplacian Resolution
    n_pairs=numel(pairs(:,1));
    c_map = linspecer(n_pairs);
    
    for i_pair = 1:n_pairs
        figure
        i_N = pairs(i_pair,1);
        i_T = pairs(i_pair,2);
        i_q = pairs(i_pair,3);

        sqrtN = sqrtN_vals(i_N);
        T = T_vals(i_T);
        q_vals = qbin{i_N,i_T};
        q=q_vals(i_q);

        gamma_cur=param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+1);

        
        t_vals=TCF_times{i_N,i_T};
        N_q=length(qbin{i_N,i_T});
        N_t=length(t_vals);
        q_indices = (i_q):N_q:N_t*N_q; 
        cf_vals=real(gmperpmperp{i_N,i_T}(q_indices));
    
        h_cf_plot=plot(t_vals,cf_vals); hold on; 
        set(h_cf_plot,'LineStyle','-','LineWidth',2,'Color','blue');
    
        xlim([0,min(1e4,30/gamma_cur)]);
        
        hXLabel = xlabel('$t$','interpreter','latex');
        hYLabel = ylabel('$C_{m\perp}(q,t)$','interpreter','latex');

        % Font
        set([gca,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize',fontsize_axis)

        pos = get(gca, 'position');
        dim = [pos(3)-.16 .8*pos(2) pos(3) pos(4)]; %top right
%     dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
%     dim = [pos(3)-.05, .1, pos(3), .3*pos(4)]; % bottom right
%     str = {dt_str,tmax_str,om_str,sig_str,tau_str};
        str = {sprintf('$N = (%d)^2$',sqrtN), sprintf('$T = %.3f$',T), sprintf('$q = %.3f$',q),'no resolution'};
        h_ann=annotation('textbox',dim,'String',str,'FitBoxToText','on',...
            'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);
        
        if(saveswitch == 1)
            figname=sprintf('%s/%s_TimeDomain_sqrtN_%d_T_%.3f_q_%.3f_NoResolution',basedir,curmodel,sqrtN,T,q);
            fprintf('Creating figure %s\n',figname)
            exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
            exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
        
        for res_factor=res_vals
            resolution=res_factor*gamma_cur;
            delete(h_cf_plot);
            h_cf_plot=plot(t_vals,cf_vals.*exp(-t_vals*resolution)); hold on; 
    
            set(h_cf_plot,'LineStyle','-','LineWidth',2,'Color','blue');

            delete(findall(gcf,'type','annotation'));
            pos = get(gca, 'position');
            dim = [pos(3)-.2 .8*pos(2) pos(3) pos(4)]; %top right
            str = {sprintf('$N = (%d)^2$',sqrtN), sprintf('$T = %.3f$',T), sprintf('$q = %.3f$',q),sprintf('$\\tau = %.3e$',1/resolution)};
            h_ann=annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);

            if(saveswitch == 1)
                figname=sprintf('%s/%s_TimeDomain_Laplacian_sqrtN_%d_T_%.3f_q_%.3f_tau_%.3e',basedir,curmodel,sqrtN,T,q,1/resolution);
                fprintf('Creating figure %s\n',figname)
                exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
            end
        end
    end


    %% 2 Behavior in frequency domain, Laplacian Resolution
    for i_pair = 1:n_pairs
        figure
        i_N = pairs(i_pair,1);
        i_T = pairs(i_pair,2);
        i_q = pairs(i_pair,3);

        sqrtN = sqrtN_vals(i_N);
        T = T_vals(i_T);
        q_vals = qbin{i_N,i_T};
        q=q_vals(i_q);

        gamma_cur=param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+1);
        omega_1_cur=param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+2);
        
        t_vals=TCF_times{i_N,i_T};
        N_q=length(qbin{i_N,i_T});
        N_t=length(t_vals);
        q_indices = (i_q):N_q:N_t*N_q; 
        cf_vals=real(gmperpmperp{i_N,i_T}(q_indices));
        
        [ft_vals,om_vals]=FT_correlation(t_vals,cf_vals, 1e6);
        h_cf_plot=plot(om_vals,ft_vals); hold on; 
        set(h_cf_plot,'LineStyle','-','LineWidth',2,'Color','blue');

        xlim([0,3*omega_1_cur]);

        hXLabel = xlabel('$\omega$','interpreter','latex');
        hYLabel = ylabel('$S_{m\perp}(q,\omega)$','interpreter','latex');

%         % Font
        set([gca,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize',fontsize_axis)

        pos = get(gca, 'position');
        dim = [pos(3)-.16 .8*pos(2) pos(3) pos(4)]; %top right
%     dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
%     dim = [pos(3)-.05, .1, pos(3), .3*pos(4)]; % bottom right
%     str = {dt_str,tmax_str,om_str,sig_str,tau_str};
        str = {sprintf('$N = (%d)^2$',sqrtN), sprintf('$T = %.3f$',T), sprintf('$q = %.3f$',q),'no resolution'};
        h_ann=annotation('textbox',dim,'String',str,'FitBoxToText','on',...
            'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);
        
        if(saveswitch == 1)
            figname=sprintf('%s/%s_OmegaDomain_sqrtN_%d_T_%.3f_q_%.3f_NoResolution',basedir,curmodel,sqrtN,T,q);
            fprintf('Creating figure %s\n',figname)
            exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
            exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
        
        for res_factor=res_vals
            resolution=res_factor*gamma_cur;
            delete(h_cf_plot);
            [ft_vals,om_vals]=FT_correlation(t_vals,cf_vals.*exp(-t_vals*resolution), 1e6);
            h_cf_plot=plot(om_vals,ft_vals); hold on; 
%             h_cf_plot=plot(t_vals,cf_vals.*exp(-t_vals*resolution)); hold on; 
    
            set(h_cf_plot,'LineStyle','-','LineWidth',2,'Color','blue');

            delete(findall(gcf,'type','annotation'));
            pos = get(gca, 'position');
            dim = [pos(3)-.2 .8*pos(2) pos(3) pos(4)]; %top right
            str = {sprintf('$N = (%d)^2$',sqrtN), sprintf('$T = %.3f$',T), sprintf('$q = %.3f$',q),sprintf('$\\tau = %.3e$',1/resolution)};
            h_ann=annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);

            if(saveswitch == 1)
                figname=sprintf('%s/%s_OmegaDomain_Laplacian_sqrtN_%d_T_%.3f_q_%.3f_tau_%.3e',basedir,curmodel,sqrtN,T,q,1/resolution);
                fprintf('Creating figure %s\n',figname)
                exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
            end
        end
    end



    %% 3 Behavior in time domain, Gaussian Resolution
    n_pairs=numel(pairs(:,1));
    c_map = linspecer(n_pairs);
    
    for i_pair = 1:n_pairs
        figure(200+i_pair)
        i_N = pairs(i_pair,1);
        i_T = pairs(i_pair,2);
        i_q = pairs(i_pair,3);

        sqrtN = sqrtN_vals(i_N);
        T = T_vals(i_T);
        q_vals = qbin{i_N,i_T};
        q=q_vals(i_q);

        gamma_cur=param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+1);

        
        t_vals=TCF_times{i_N,i_T};
        N_q=length(qbin{i_N,i_T});
        N_t=length(t_vals);
        q_indices = (i_q):N_q:N_t*N_q; 
        cf_vals=real(gmperpmperp{i_N,i_T}(q_indices));
    
        h_cf_plot=plot(t_vals,cf_vals); hold on; 
        set(h_cf_plot,'LineStyle','-','LineWidth',2,'Color','blue');
    
        xlim([0,min(1e4,30/gamma_cur)]);
        
        hXLabel = xlabel('$t$','interpreter','latex');
        hYLabel = ylabel('$C_{m\perp}(q,t)$','interpreter','latex');

        % Font
        set([gca,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize',fontsize_axis)

        pos = get(gca, 'position');
        dim = [pos(3)-.16 .8*pos(2) pos(3) pos(4)]; %top right
%     dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
%     dim = [pos(3)-.05, .1, pos(3), .3*pos(4)]; % bottom right
%     str = {dt_str,tmax_str,om_str,sig_str,tau_str};
        str = {sprintf('$N = (%d)^2$',sqrtN), sprintf('$T = %.3f$',T), sprintf('$q = %.3f$',q),'no resolution'};
        h_ann=annotation('textbox',dim,'String',str,'FitBoxToText','on',...
            'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);
        
        
        for res_factor=res_vals
            resolution=res_factor*gamma_cur;
            delete(h_cf_plot);
            h_cf_plot=plot(t_vals,cf_vals.*exp(-.5*t_vals.^2*resolution^2)); hold on; 
    
            set(h_cf_plot,'LineStyle','-','LineWidth',2,'Color','blue');

            delete(findall(gcf,'type','annotation'));
            pos = get(gca, 'position');
            dim = [pos(3)-.2 .8*pos(2) pos(3) pos(4)]; %top right
            str = {sprintf('$N = (%d)^2$',sqrtN), sprintf('$T = %.3f$',T), sprintf('$q = %.3f$',q),sprintf('$\\sigma = %.3e$',1/resolution)};
            h_ann=annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);

            if(saveswitch == 1)
                figname=sprintf('%s/%s_TimeDomain_Gaussian_sqrtN_%d_T_%.3f_q_%.3f_sigma_%.3e',basedir,curmodel,sqrtN,T,q,1/resolution);
                fprintf('Creating figure %s\n',figname)
                exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
            end
        end
    end


    %% 4 Behavior in frequency domain, Gaussian Resolution
    for i_pair = 1:n_pairs
        figure(300+i_pair)
        i_N = pairs(i_pair,1);
        i_T = pairs(i_pair,2);
        i_q = pairs(i_pair,3);

        sqrtN = sqrtN_vals(i_N);
        T = T_vals(i_T);
        q_vals = qbin{i_N,i_T};
        q=q_vals(i_q);

        gamma_cur=param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+1);
        omega_1_cur=param_TCFSpin_omega_1_q_DO{i_N,i_T}(2*i_q+2);
        
        t_vals=TCF_times{i_N,i_T};
        N_q=length(qbin{i_N,i_T});
        N_t=length(t_vals);
        q_indices = (i_q):N_q:N_t*N_q; 
        cf_vals=real(gmperpmperp{i_N,i_T}(q_indices));
        
        [ft_vals,om_vals]=FT_correlation(t_vals,cf_vals, 1e6);
        h_cf_plot=plot(om_vals,ft_vals); hold on; 
        set(h_cf_plot,'LineStyle','-','LineWidth',2,'Color','blue');

        xlim([0,3*omega_1_cur]);

        hXLabel = xlabel('$\omega$','interpreter','latex');
        hYLabel = ylabel('$S_{m\perp}(q,\omega)$','interpreter','latex');

%         % Font
        set([gca,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize',fontsize_axis)

        pos = get(gca, 'position');
        dim = [pos(3)-.16 .8*pos(2) pos(3) pos(4)]; %top right
%     dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
%     dim = [pos(3)-.05, .1, pos(3), .3*pos(4)]; % bottom right
%     str = {dt_str,tmax_str,om_str,sig_str,tau_str};
        str = {sprintf('$N = (%d)^2$',sqrtN), sprintf('$T = %.3f$',T), sprintf('$q = %.3f$',q),'no resolution'};
        h_ann=annotation('textbox',dim,'String',str,'FitBoxToText','on',...
            'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);
      
        for res_factor=res_vals
            resolution=res_factor*gamma_cur;
            delete(h_cf_plot);
            [ft_vals,om_vals]=FT_correlation(t_vals,cf_vals.*exp(-.5*t_vals.^2*resolution^2), 1e6);
            h_cf_plot=plot(om_vals,ft_vals); hold on; 
%             h_cf_plot=plot(t_vals,cf_vals.*exp(-t_vals*resolution)); hold on; 
    
            set(h_cf_plot,'LineStyle','-','LineWidth',2,'Color','blue');

            delete(findall(gcf,'type','annotation'));
            pos = get(gca, 'position');
            dim = [pos(3)-.2 .8*pos(2) pos(3) pos(4)]; %top right
            str = {sprintf('$N = (%d)^2$',sqrtN), sprintf('$T = %.3f$',T), sprintf('$q = %.3f$',q),sprintf('$\\sigma = %.3e$',1/resolution)};
            h_ann=annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);

            if(saveswitch == 1)
                figname=sprintf('%s/%s_OmegaDomain_Gaussian_sqrtN_%d_T_%.3f_q_%.3f_sigma_%.3e',basedir,curmodel,sqrtN,T,q,1/resolution);
                fprintf('Creating figure %s\n',figname)
                exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
            end
        end
    end
end
    