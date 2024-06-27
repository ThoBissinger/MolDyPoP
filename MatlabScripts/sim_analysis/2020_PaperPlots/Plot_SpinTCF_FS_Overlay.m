% clear
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/TCF/mperp_TCF_Overlay',fig_base);

for i_model = [3]

    % figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=matfile(mxydata_AdjustedTime_name);
        fitdata=matfile(mxyfit_TCF_q_name);
        
        i_T = 4;
        % pairs of i_N, i_T, i_q
%         pairs=[2,10,1; 3,10,3; 4,10,9; 5,10,23];
        pairs=[3,i_T,1; 4,i_T,3; 5,i_T,9];
        
        t_max = 3e3;
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
        
%         data=matfile(fmxydata_AdjustedTime_name);
%         fitdata=matfile(fmxyfit_TCF_q_name);

        i_T = 4; t_max = 5.5e3;
%         i_T=10; t_max = 3e3;
        pairs=[3,i_T,1; 4,i_T,3; 5,i_T,9];

        

        L_vals=[9.25,18.5,37,74,148];
        runmax=500;
    end
    fontsize_axis=15;
    fontsize_annotation=15;

    
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
    figure
    n_pairs=numel(pairs(:,1));
    c_map = linspecer(n_pairs);
    hold on;
    for i_pair = 1:n_pairs
        i_N = pairs(i_pair,1);
%         i_T = pairs(i_pair,2);
        i_q = pairs(i_pair,3);

        sqrtN = sqrtN_vals(i_N);
        T = T_vals(i_T);
        q_vals = qbin{i_N,i_T};
        q=q_vals(i_q);
        dispname=sprintf('$N = (%d)^2$',sqrtN);
        
        [h1,h2] = CreateSingle_TCF_Plot(TCF_times,gmperpmperp,qbin,i_N,i_T,i_q);
        delete(h2);
        set(h1,'YData',h1.YData/h1.YData(1));
        set(h1,'LineStyle','-','LineWidth',2,'Color',c_map(i_pair,:),...
            'DisplayName',dispname);
    end
    xlim([0,t_max]);
    hXLabel = xlabel('$t$','interpreter','latex');
    hYLabel = ylabel('$C_{m\perp}(q,t)/\chi_{m\perp}(q)$','interpreter','latex');
    hlegend = legend('Location', 'NorthEast', ...
        'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);
    set([gca,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize',fontsize_axis);


    pos = get(gca, 'position');
%     dim = [pos(1) pos(2) 1.5*pos(3) .5*pos(4)]; %bottom right
%     dim = [ .6 .15 .9 .5];
    dim = [ pos(3)-.75*pos(1), pos(2), .9, .15]; %bottom right
    str = {sprintf('$T = %.3f$',T), sprintf('$q = %.3f$',q)};
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation);
    if(saveswitch == 1)
        figname=sprintf('%s/%s_SpinTCF_Overlay_T_%.3f_q_%.3f',basedir,curmodel,T,q);
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end
    