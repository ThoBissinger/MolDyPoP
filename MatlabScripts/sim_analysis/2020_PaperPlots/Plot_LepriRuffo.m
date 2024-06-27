%% 0 Initialize
run initialization_script;
saveswitch=0;
basedir=sprintf('%s/plots/LepriRuffo',fig_base);
fontsize_axis=15;
fontsize_annotation=15;

section_select = 1;
    
for i_model = [1]
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";

        data=matfile(mxydata_LepriRuffo_name);
        data_LongTime=matfile(mxydata_LepriRuffo_name);
        fitdata=matfile(mxyfit_name);
        fitdata_FS=matfile(mxyfit_FS_name);
        
        eta_dataset_id='dynamics_mxy_LinearTime';
%         T_select=[8:16];
        T_select=[1:19];
        T_select=[1:12];
        L_vals=[9.25,18.5,37,74,148];
        
        y_offsets_fit=[.97,.9,.87];

        t_max_long = 8;
        t_max_short = 1;
        
    elseif (i_model == 2)
        curmodel="xy";
        curtitle="XY model";

        data=matfile(xydata_LepriRuffo_name);
        data_LongTime=matfile(xydata_LepriRuffo_name);
        fitdata=matfile(xyfit_name);
        fitdata_FS=matfile(xyfit_FS_name);
                
        T_select=[8:16];
        L_vals=[16, 32, 64, 128, 256];

        eta_best=fitdata_FS.eta_vals;
        y_offsets_fit=1 - eta_best;
        % eta_vals = [0.0069, 0.0230, 0.0383];
        
        t_max_long = 6;
        t_max_short = .3;

    elseif (i_model == 3)
        curmodel="fmxy";
        curtitle="FMXY model";

        data=matfile(fmxydata_LepriRuffo_name);
        data_LongTime=matfile(fmxydata_LepriRuffo_name);
        fitdata=matfile(fmxyfit_name);
        fitdata_FS=matfile(fmxyfit_FS_name);

        eta_dataset_id='dynamics_fmxy_LinearTime';
        
        T_select=[8:16];
        T_select=[17:20];
        T_select=1:11;
        L_vals=[9.25,18.5,37,74,148];
        
        y_offsets_fit=[.97,.9,.87];

        t_max_long = 10;
        t_max_short = 2;
    end

    T_vals=data.('T_vals');
    sqrtN_vals=data.('sqrtN_vals');
    i_N = numel(sqrtN_vals);
    
    ACF_Spin = data.('ACF_Spin');
    ACF_Spin_LinearTime = data_LongTime.('ACF_Spin');
    absM_av = data_LongTime.('absM_av');
    
    averaging_times = data.('averaging_times');
    averaging_times_LinearTime = data_LongTime.('averaging_times');
    t_vals = cell2mat(averaging_times(end,1));

    eta_best=fitdata_FS.eta_vals;
    
    gui_locations={'northwest', 'north', 'northeast','west','center','east','southwest','south','southeast'};
    %     [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% section 1 Create ACF plot for long time (Linear Time sampling)
    c_map = linspecer(numel(sqrtN_vals));
    if (ismember(1,section_select))
        for i = 1:numel(T_select)
            fig_ind = 1000 * i_model + i;
            curfig=figure;
            movegui(curfig,gui_locations{min(i,numel(gui_locations))})
            hold on
            i_T = T_select(i);
            T = T_vals(i_T);
            i_T_eta = find_T_index(eta_dataset_id,T);
            eta_cur=eta_best(i_T_eta);
        %         color_ind = floor(T / T_max * T_grating);
            for i_N = 1:numel(sqrtN_vals)
    %             curACF=cell2mat(ACF_Spin(i_N,i_T));
                curACF=ACF_Spin_LinearTime{i_N,i_T};
                t_vals=averaging_times_LinearTime{i_N,i_T};
    %             curabsM=cell2mat(absM_av(i_N,i_T));
                dispname = sprintf('$N = %i^2$',sqrtN_vals(i_N));
    %             semilogx(t_vals/L_vals(i_N), curACF, ...
    %                 'LineStyle', '-', 'LineWidth', 2, ...
    %                 'Displayname',dispname,...
    %                 'Color',c_map(i_N,:));
    %             figure(i)

%                 hACF_line{i} = plot(t_vals/L_vals(i_N), L_vals(i_N) ^ eta_cur * curACF);
                hACF_line{i} = plot(t_vals/L_vals(i_N), curACF / absM_av{i_N,i_T}^2);
                set(hACF_line{i}, ...
                    'LineStyle', '-', 'LineWidth', 2, ...
                    'Displayname',dispname,...
                    'Color',c_map(i_N,:));
            end

            % Legend, axes etc
            hLegend{i} = legend('Location', 'NorthEast','interpreter','latex','FontSize', fontsize_annotation,...
            'NumColumns',1);

            hXLabel{i} = xlabel('$t / L$','interpreter','latex');
    %         hYLabel{i} = ylabel('$L^{\eta} \langle \mathbf{s}_i(t) \cdot \mathbf{s}_i(0)\rangle$','interpreter','latex');
            hYLabel{i} = ylabel('$L^{\eta} C_m^{\textrm{inc}}(t)$','interpreter','latex');
            hTitle{i} = title(curtitle);

            cur_gca{i} = gca;
            xlim([0 t_max_long]);
            yl = ylim;
            ymax = yl(2);
            % Font
            set([cur_gca{i},hXLabel{i}, hYLabel{i}], 'FontName', 'cmr12', 'FontSize', fontsize_axis)



            % Adjust axes properties
            set(cur_gca{i}, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
                'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
                'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
                'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')


            % Text

            text(0.05,.95,sprintf('$T = %.3f$',T),...
                'Units','normalized',...
                'VerticalAlignment','top', 'HorizontalAlignment','left',...
                'interpreter','latex',...
                'Color','Red','FontSize', 10);

            text(0.05,.90,sprintf('$\\eta = %.4f$',eta_cur),...
                'Units','normalized',...
                'VerticalAlignment','top', 'HorizontalAlignment','left',...
                'interpreter','latex',...
                'Color','Red','FontSize', 10);



            if(saveswitch == 1)
                figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/LepriRuffo/%s_LepriRuffo_T_%.3f',curmodel,T);
                fprintf('Creating figure %s\n',figname)
                exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
            end

        end
    end
    
    
    
    %% 2 Create ACF plot for short time
    c_map = linspecer(numel(sqrtN_vals));
    if (ismember(2,section_select))
        for i = 1:numel(T_select)
            fig_ind = 1000 * i_model + i + 100;
            curfig=figure;
            movegui(curfig,gui_locations{min(i,numel(gui_locations))})
            hold on
            i_T = T_select(i);
            T = T_vals(i_T);
            i_T_eta = find_T_index(eta_dataset_id,T);
            eta_cur=eta_best(i_T_eta);
        %         color_ind = floor(T / T_max * T_grating);
            for i_N = 1:numel(sqrtN_vals)
    %             curACF=cell2mat(ACF_Spin(i_N,i_T));
                curACF=ACF_Spin{i_N,i_T};
                t_vals=averaging_times{i_N,i_T};
    %             curabsM=cell2mat(absM_av(i_N,i_T));
                dispname = sprintf('$N = %i^2$',sqrtN_vals(i_N));
    %             semilogx(t_vals/L_vals(i_N), curACF, ...
    %                 'LineStyle', '-', 'LineWidth', 2, ...
    %                 'Displayname',dispname,...
    %                 'Color',c_map(i_N,:));
    %             figure(i)

                hACF_line{i} = plot(t_vals/L_vals(i_N), L_vals(i_N) ^ eta_cur * curACF);
    %             hACF_line{i} = plot(t_vals, curACF);
                set(hACF_line{i}, ...
                    'LineStyle', '-', 'LineWidth', 2, ...
                    'Displayname',dispname,...
                    'Color',c_map(i_N,:));
            end

            % Legend, axes etc
            hLegend{i} = legend('Location', 'NorthEast','interpreter','latex','FontSize', 10,...
            'NumColumns',1);

            hXLabel{i} = xlabel('$t / L$','interpreter','latex');
            hYLabel{i} = ylabel('$L^{\eta} C_m^{\textrm{inc}}(t)$','interpreter','latex');
            hTitle{i} = title(curtitle);

            cur_gca{i} = gca;

            yl = ylim;
            ymax = yl(2);
            % Font
            set(cur_gca{i}, 'FontName', 'cmr12')
            set([hXLabel{i}, hYLabel{i}], 'FontName', 'cmr12')
            set(cur_gca{i}, 'FontSize', 8)
            set([hXLabel{i}, hYLabel{i}], 'FontSize', 12)
            set(hTitle{i}, 'FontSize', 12, 'FontWeight' , 'bold')



            % Adjust axes properties
            set(cur_gca{i}, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
                'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
                'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
                'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')


            % Text

            text(0.05,.1,sprintf('$T = %.2f$',T),...
                'Units','normalized',...
                'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
                'interpreter','latex',...
                'Color','Red','FontSize', 10);

            text(0.05,.05,sprintf('$\\eta = %.4f$',eta_cur),...
                'Units','normalized',...
                'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
                'interpreter','latex',...
                'Color','Red','FontSize', 10);




%             shorttimefig=figure;
%             movegui(figure,gui_locations{min(i,numel(gui_locations))})
    %         ax_shorttimefig=copyobj(ax_longtimefig,shorttimefig);
            legend show;
            legend('Location', 'NorthEast','interpreter','latex','FontSize', 10,...
            'NumColumns',1);
        set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
                'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
                'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
                'LineWidth', .5,'Xscale', 'log','Yscale', 'log')
    %         movegui(figure(10 * i_model + i),'Position',[0 1], 'units','normalized')
            xx=linspace(.1,.8);
    %         dispname='$\sim \eta \ln(t/L)$';
    %         plot(xx,y_offsets_fit(i_T)-eta_best(i_T)*log(xx),'--','DisplayName',dispname,'Color','black');
            dispname='$\sim (t/L)^{-\eta}$';
            plot(xx,i_T_eta^2*.001+xx.^-eta_best(i_T_eta),'--','DisplayName',dispname,...
                'Color','black','LineWidth',2);
            xlim([1e-2 t_max_short])
    %         ylim([.9 1.1])
            if(saveswitch == 1)
                figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/LepriRuffo/%s_LepriRuffo_ShortTime_T_%.3f',curmodel,T);
                fprintf('Creating figure %s\n',figname)
                exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
            end
        end
    end
    
    %% 3 Create ACF plot for short time with better scaling
    c_map = linspecer(numel(sqrtN_vals));
    if (ismember(3,section_select))
        for i = 1:numel(T_select)
            fig_ind = 1000 * i_model + i + 200;
            curfig=figure;
            movegui(curfig,gui_locations{min(i,numel(gui_locations))})
            hold on
            i_T = T_select(i);
            T = T_vals(i_T);
            i_T_eta = find_T_index(eta_dataset_id,T);
            eta_cur=eta_best(i_T_eta);
        %         color_ind = floor(T / T_max * T_grating);
            for i_N = 1:numel(sqrtN_vals)
    %             curACF=cell2mat(ACF_Spin(i_N,i_T));
                curACF=ACF_Spin{i_N,i_T};
                t_vals=averaging_times{i_N,i_T};
    %             curabsM=cell2mat(absM_av(i_N,i_T));
                dispname = sprintf('$N = %i^2$',sqrtN_vals(i_N));
    %             semilogx(t_vals/L_vals(i_N), curACF, ...
    %                 'LineStyle', '-', 'LineWidth', 2, ...
    %                 'Displayname',dispname,...
    %                 'Color',c_map(i_N,:));
    %             figure(i)
                t_exponent = 0;
                C_exponent = 0;
                hACF_line{i} = plot(t_vals/L_vals(i_N)^t_exponent, L_vals(i_N) ^ C_exponent * curACF);
    %             hACF_line{i} = plot(t_vals, curACF);
                set(hACF_line{i}, ...
                    'LineStyle', '-', 'LineWidth', 2, ...
                    'Displayname',dispname,...
                    'Color',c_map(i_N,:));
            end

            % Legend, axes etc
            hLegend{i} = legend('Location', 'NorthWest','interpreter','latex','FontSize', 10,...
            'NumColumns',1);

            hXLabel{i} = xlabel('$t$','interpreter','latex');
            hYLabel{i} = ylabel('$C_m^{\textrm{inc}}(t)$','interpreter','latex');
            hTitle{i} = title(curtitle);

            cur_gca{i} = gca;

            yl = ylim;
            ymax = yl(2);
            % Font
            set(cur_gca{i}, 'FontName', 'cmr12')
            set([hXLabel{i}, hYLabel{i}], 'FontName', 'cmr12')
            set(cur_gca{i}, 'FontSize', 8)
            set([hXLabel{i}, hYLabel{i}], 'FontSize', 12)
            set(hTitle{i}, 'FontSize', 12, 'FontWeight' , 'bold')



            % Adjust axes properties
            set(cur_gca{i}, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
                'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
                'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
                'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')


            % Text

            text(0.05,.1,sprintf('$T = %.2f$',T),...
                'Units','normalized',...
                'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
                'interpreter','latex',...
                'Color','Red','FontSize', 10);

            text(0.05,.05,sprintf('$\\eta = %.4f$',eta_cur),...
                'Units','normalized',...
                'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
                'interpreter','latex',...
                'Color','Red','FontSize', 10);




%             curfig=figure;
%             movegui(curfig,gui_locations{min(i,numel(gui_locations))})
    %         ax_shorttimefig=copyobj(ax_longtimefig,shorttimefig);
            legend show;
            legend('Location', 'NorthEast','interpreter','latex','FontSize', 10,...
            'NumColumns',1);
        set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
                'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
                'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
                'LineWidth', .5,'Xscale', 'log','Yscale', 'log')
    %         movegui(figure(10 * i_model + i),'Position',[0 1], 'units','normalized')
            xx=linspace(1,sqrtN_vals(end));
    %         dispname='$\sim \eta \ln(t/L)$';
    %         plot(xx,y_offsets_fit(i_T)-eta_best(i_T)*log(xx),'--','DisplayName',dispname,'Color','black');
            dispname='$\sim (t)^{-\eta}$';
            plot(xx,i_T_eta^2*.001+.05+xx.^-eta_best(i_T_eta),'--','DisplayName',dispname,'Color','black');
            xlim([0 1.5*sqrtN_vals(end)])
    %         ylim([.9 1.1])
            if(saveswitch == 1)
                figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/LepriRuffo/%s_LepriRuffo_ShortTime_Bare_T_%.3f',curmodel,T);
                fprintf('Creating figure %s\n',figname)
                exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
                exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
            end
        end
    end
end
   