run initialization_script;
saveswitch=0;
basedir=sprintf('%s/plots/TCF/mperp_TCF_WithFit',fig_base);
fontsize_axis=15;
fontsize_annotation=15;
for i_model = [1]

    % figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=matfile(mxydata_name);
        fitdata=matfile(mxyfit_TCF_q_name);
        
        T_select=[9:2:13,14:21];
        T_select=1:24;
        
        L_vals=[9.25,18.5,37,74,148];
        
        q_select=1:2:6;
        
    elseif (i_model == 2)
        curmodel="xy";
        curtitle="SXY model";
        
        data=matfile(xydata_name);
        fitdata=matfile(xyfit_TCF_q_name);
        
        T_select=[17:3:38];
        
        L_vals=[16,32,64,128,256];
        
        q_select=1:6;
    elseif (i_model == 3)
        curmodel="fmxy";
        curtitle="FMXY model";
        
        data=matfile(fmxydata_name);
        fitdata=matfile(fmxyfit_TCF_q_name);
        
        T_select=[9:2:13,14:21];
        T_select=1:24;
        
        L_vals=[9.25,18.5,37,74,148];
        
        q_select=1:2:6;
        
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
    %% 1 Without fit
    N_sqrtN = numel(sqrtN_vals);
    c_map = linspecer(numel(q_select));
    
    for img_index = 1:length(T_select)
        figure(img_index)
        i_T = T_select(img_index);
        T = T_vals(i_T);
        
        gmperpmperp_cur = gmperpmperp{i_N,i_T};
        
        ymax=0;
        ymin=0;
        for i_q = 1:numel(q_select)
            index_q = q_select(i_q);
            q = q_vals(index_q);
            dispname=sprintf('q = %.3f', q);
    %     coeffs = cell2mat(param_TCFSpin_omega_1_q_DO(i_N,:)');
        

            
            t_vals=TCF_times{i_N,i_T};
            TCF_Spin_cur = real(gmperpmperp_cur(index_q:length(q_vals):end));

            plot(t_vals,TCF_Spin_cur/TCF_Spin_cur(1),...
                'Color',c_map(i_q,:),'DisplayName',dispname,...
                'LineStyle', '-', ...
                'LineWidth',1.5);
%                 'Marker', '^', 'MarkerSize', 4, ...
%                 'LineWidth',1.5); 
            hold on;
            
            ymax=max(ymax,max(abs(TCF_Spin_cur)));
            ymin=min(ymin,min(real(TCF_Spin_cur)));
    %         legend show
    %         legend('Location','best')
        end
%         ymax=1.2;
%         ymin=-1.2;
        h_axis = gca;
        set(h_axis,'xscale','log');

        hXLabel = xlabel('$t$','interpreter','latex');
        hYLabel = ylabel('$C_{m\perp}(q,t)/\chi_{m\perp}$','interpreter','latex');



%         ylim([min(1.1*ymin,-.4*ymax),1.3*ymax]);
        ylim([-1.2,2]);
        xlim([t_vals(1),2e3]);

        % Legend, axes etc
        hLegend = legend('Location', 'NorthWest','interpreter','latex','FontSize', fontsize_annotation,...
            'NumColumns',1);

        hTitle = title(curtitle);

        % Font
%         set(h_axis, 'FontName', 'cmr12')
%         set([hXLabel, hYLabel], 'FontName', 'cmr12')
%         set(h_axis, 'FontSize', 12)
        set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12', 'FontSize', fontsize_axis)
%         set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


        % Adjust axes properties
        set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
            'LineWidth', .5,'Xscale', 'log','Yscale', 'lin')
        
        pos = get(gca, 'position');
%     dim = [pos(1) pos(2) 1.5*pos(3) .5*pos(4)]; %bottom right
        dim = [pos(3)-.1 .8*pos(2) pos(3) pos(4)]; %top right
%   dim = [ pos(3)-.75*pos(1), pos(2), .9, .15]; %bottom right
        str = {sprintf('$N = %d^2$',sqrtN_vals(i_N)),sprintf('$T = %.3f$',T)};
        annotation('textbox',dim,'String',str,'FitBoxToText','on',...
            'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation,...
            'BackgroundColor','white');
%         text(0.95,.95,sprintf('$N = %d^2$',sqrtN_vals(i_N)),...
%             'Units','normalized',...
%             'VerticalAlignment','top', 'HorizontalAlignment','right',...
%             'interpreter','latex',...
%             'Color','Blue','FontSize', 10);
% 
%         text(0.95,.90,sprintf('$T = %.3f$',T),...
%             'Units','normalized',...
%             'VerticalAlignment','top', 'HorizontalAlignment','right',...
%             'interpreter','latex',...
%             'Color','Blue','FontSize', 10);

        figname=sprintf('%s/%s_SpinTCF_sqrtN_%d_T_%.3f',basedir,curmodel,sqrtN_vals(i_N),T);
        fprintf('Creating figure %s\n',figname)
        if(saveswitch == 1)
           exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
        
        xlim([t_vals(1),t_vals(end)]);
        set(h_axis, 'Xscale', 'linear');
        figname=sprintf('%s/%s_SpinTCF_Linear_sqrtN_%d_T_%.3f',basedir,curmodel,sqrtN_vals(i_N),T);
        fprintf('Creating figure %s\n',figname)
        if(saveswitch == 1)
           exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
    end
    
    %% 2 With fits
    for img_index = 1:length(T_select)
        figure(img_index+100)
        i_T = T_select(img_index);
        T = T_vals(i_T);
        
        gmperpmperp_cur = gmperpmperp{i_N,i_T};
        
        curcoeff = param_TCFSpin_omega_1_q_DO{i_N,i_T};
        
        ymax=0;
        ymin=0;
        for i_q = 1:numel(q_select)
            index_q = q_select(i_q);
            q = q_vals(index_q);
            dispname=sprintf('q = %.3f', q);
    %     coeffs = cell2mat(param_TCFSpin_omega_1_q_DO(i_N,:)');
        

            
            t_vals=TCF_times{i_N,i_T};
            TCF_Spin_cur = gmperpmperp_cur(index_q:length(q_vals):end);

            plot(t_vals,real(TCF_Spin_cur/TCF_Spin_cur(1)),...
                'Color',c_map(i_q,:),...
                'LineStyle', 'none', 'HandleVisibility', 'off',... % 'IconDisplayStyle', 'off', ...
                'Marker', 'o', 'MarkerSize', 4, 'LineWidth',1.8,...
                'MarkerIndices',unique(floor(logspace(0,3,300))));
            hold on;
                        
            c=curcoeff(2*(index_q)+(1:2));
%             c=fit_DampedOscillator_RealSpace(t_vals,TCF_Spin_cur,6,1,'omega_1');
            plot(t_vals,fitfunc_DO(t_vals,1,c),...
                'Color',c_map(i_q,:),'DisplayName',dispname,...
                'LineStyle', '-', ...
                'LineWidth',1.8);
            
%             ymax=max(ymax,max(abs(TCF_Spin_cur)));
%             ymin=min(ymin,min(real(TCF_Spin_cur)));
    %         legend show
    %         legend('Location','best')
        end

        h_axis = gca;
        set(h_axis,'xscale','log');

        hXLabel = xlabel('$t$','interpreter','latex');
        hYLabel = ylabel('$C_{m\perp}(q,t)/\chi_{m\perp}$','interpreter','latex');



        ylim([-1.2,2]);
        xlim([t_vals(1),t_vals(end)]);

        % Legend, axes etc
        hLegend = legend('Location', 'NorthWest','interpreter','latex','FontSize', fontsize_annotation,...
            'NumColumns',1);

        hTitle = title(curtitle);

        % Font
        set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12', 'FontSize', fontsize_axis)


        % Adjust axes properties
        set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
            'LineWidth', .5,'Xscale', 'log','Yscale', 'lin')

        pos = get(gca, 'position');
%     dim = [pos(1) pos(2) 1.5*pos(3) .5*pos(4)]; %bottom right
        dim = [pos(3)-.1 .8*pos(2) pos(3) pos(4)]; %top right
%       dim = [ pos(3)-.75*pos(1), pos(2), .9, .15]; %bottom right
        str = {sprintf('$N = %d^2$',sqrtN_vals(i_N)),sprintf('$T = %.3f$',T)};
        annotation('textbox',dim,'String',str,'FitBoxToText','on',...
            'interpreter','latex','FontName', 'cmr12','FontSize',fontsize_annotation,...
            'BackgroundColor','white');
%         text(0.95,.95,sprintf('$N = %d^2$',sqrtN_vals(i_N)),...
%             'Units','normalized',...
%             'VerticalAlignment','top', 'HorizontalAlignment','right',...
%             'interpreter','latex',...
%             'Color','Blue','FontSize', 10);
% 
%         text(0.95,.90,sprintf('$T = %.3f$',T),...
%             'Units','normalized',...
%             'VerticalAlignment','top', 'HorizontalAlignment','right',...
%             'interpreter','latex',...
%             'Color','Blue','FontSize', 10);

        figname=sprintf('%s/%s_SpinTCF_WithFit_sqrtN_%d_T_%.3f',basedir,curmodel,sqrtN_vals(i_N),T);
        fprintf('Creating figure %s\n',figname)
        if(saveswitch == 1)
           exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
        
        set(h_axis, 'Xscale', 'linear');
        figname=sprintf('%s/%s_SpinTCF_WithFit_Linear_sqrtN_%d_T_%.3f',basedir,curmodel,sqrtN_vals(i_N),T);
        fprintf('Creating figure %s\n',figname)
        if(saveswitch == 1)
           exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
    end
end
    