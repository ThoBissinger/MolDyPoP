run initialization_script;
saveswitch=1;

% mxydata=load('mxy/rho_3.00_integ.mat');
% mxydata=load('mxy/rho_3.00_dynamics_LinearTime.mat'); mxyfit=load('mxy/rho_3.00_TCF_q_fit_LinearTime.mat');
% xydata=load('xy/lf0_qreduced.mat'); xyfit=load('xy/lf0_CritExpFit.mat');


for i_model = [1]

    % figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=matfile(mxydata_name);
        fitdata=matfile(mxyfit_TCF_q_name);
        
        % pairs of i_N, i_T, i_q
        pairs=[4,10,3; 3,10,1; 5,10,6];
        
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
        
        pairs=[4,10,3; 3,10,1; 5,10,6];

        L_vals=[9.25,18.5,37,74,148];
        runmax=500;
    end
    fontsize_axis=15;
    fontsize_annotation=15;

    basedir=sprintf('%s/plots/TCF/mperp_TCF_Anomalous',fig_base);
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
    %% 1 Long-time behavior
    n_pairs=numel(pairs(:,1));
    c_map = linspecer(n_pairs);
    
    for i_pair = 1:n_pairs
        figure(i_pair)
        i_N = pairs(i_pair,1);
        i_T = pairs(i_pair,2);
        i_q = pairs(i_pair,3);

        sqrtN = sqrtN_vals(i_N);
        T = T_vals(i_T);
        q_vals = qbin{i_N,i_T};
        q=q_vals(i_q);

        
        [h1,h2] = CreateSingle_TCF_Plot(TCF_times,gmperpmperp,qbin,i_N,i_T,i_q);
        delete(h2);
        y_errthresh=max(h1.YData)/sqrt(runmax);
        set(h1,'LineStyle','-','LineWidth',2,'Color','blue');
        yline(y_errthresh,'LineStyle',':','LineWidth',1.5,'Color','red');
        yline(-y_errthresh,'LineStyle',':','LineWidth',1.5,'Color','red');
        hXLabel = xlabel('$t$','interpreter','latex');
        hYLabel = ylabel('$C_{m\perp}(q,t)$','interpreter','latex');

%         % Font
        set([gca,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize',fontsize_axis)

        pos = get(gca, 'position');
        dim = [pos(3)-.13 .8*pos(2) pos(3) pos(4)]; %top right
%     dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
%     dim = [pos(3)-.05, .1, pos(3), .3*pos(4)]; % bottom right
%     str = {dt_str,tmax_str,om_str,sig_str,tau_str};
        str = {sprintf('$N = (%d)^2$',sqrtN), sprintf('$T = %.3f$',T), sprintf('$q = %.3f$',q)};
        annotation('textbox',dim,'String',str,'FitBoxToText','on',...
            'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);
        
         text(.96*1e4,-2*y_errthresh,'$$\frac{\chi_{m\perp}}{\sqrt{\#runs}}$$',...
            'VerticalAlignment','top','HorizontalAlignment','right',...
            'interpreter','latex',...
            'Color','red','FontSize', 12);
        if(saveswitch == 1)
            figname=sprintf('%s/%s_SpinTCF_Anomalous_sqrtN_%d_T_%.3f_q_%.3f',basedir,curmodel,sqrtN,T,q);
            fprintf('Creating figure %s\n',figname)
            exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
            exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        

            xlim([4e3,7e3]);
            ylim(4*y_errthresh*[-1,1]);
            figname=sprintf('%s/%s_SpinTCF_Anomalous_sqrtN_%d_T_%.3f_q_%.3f_zoomed',basedir,curmodel,sqrtN,T,q);
            fprintf('Creating figure %s\n',figname)
            exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
            exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
    end
end
    