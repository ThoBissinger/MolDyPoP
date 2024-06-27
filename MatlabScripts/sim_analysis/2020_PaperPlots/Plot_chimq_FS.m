clear all
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;

mxydata=load('mxy/rho_3.00_dynamics_LinearTime.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat'); mxyfit_FS=load('mxy/rho_3.00_DataCollapse_SCF.mat');
xydata=load('xy/xy_dynamics_LinearTime.mat'); xyfit=load('xy/xy_CritExpFit.mat'); xyfit_FS=load('xy/xy_DataCollapse_SCF.mat');
% xydata=load('xy/lf0_eq.mat'); xyfit=load('xy/lf0_CritExpFit.mat');


for i_model = 1:1

    figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=mxydata;
        fitdata=mxyfit;
        
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        absM_av=cell2mat(data.('absM_av'));
        SCF_Spin_av = data.('SCF_Spin_av');
        rbin = data.('rbin');
        qbin = data.('qbin');
        chimparq_av=data.('chimparq_av');
        chimperpq_av=data.('chimperpq_av');
        chimxq_av=data.('chimxq_av');
        chimyq_av=data.('chimyq_av');
        chiteq_av=data.('chiteq_av');
        chiwq_av=data.('chiwq_av');
        
        
        
%         plot_type={"mpar", "mperp", "m", "w", "te"};
        plot_type={"mpar", "mperp", "m", "te"};
%         plot_type={"mperp"};
        
        expfit=fitdata.('param_SCFSpin_Exp');
        powfit=fitdata.('param_SCFSpin_Pow');
        eta_vals_fit=fitdata.('eta_vals');
        eta_vals_FS = mxyfit_FS.('eta_vals');
        xi_vals=fitdata.('xi_vals');
        
        i_N = numel(sqrtN_vals);
        r_vals = rbin{i_N,1};
        
%         T_select=1:numel(T_vals);
%         T_select=[2,3,4,5,7];
%         T_select=[9:13,14:2:21];
        T_select=[1:19];
        
        L_vals=[9.25,18.5,37,74,148];
        r_min = 2;
        r_max = 60;
        
        T_max = .4;
        FSplot_min = .13;
        FSplot_max = .23;
        
    else
        curmodel="xy";
        curtitle="SXY model";
        
        data=xydata;
        fitdata=xyfit;
        
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        absM_av=cell2mat(data.('absM_av'));
        SCF_Spin_av = data.('SCF_Spin_av');
        rbin = data.('rbin');
        qbin = data.('qbin');
        chimparq_av=data.('chimparq_av');
        chimperpq_av=data.('chimperpq_av');
        chimxq_av=data.('chimxq_av');
        chimyq_av=data.('chimyq_av');
        chiteq_av=data.('chiteq_av');
        chiwq_av=data.('chiwq_av');
        
        
        
%         plot_type={"mpar", "mperp", "m", "w", "te"};
        plot_type={"mpar", "mperp", "m", "te"};
%         plot_type={"mperp"};
        
        expfit=fitdata.('param_SCFSpin_Exp');
        powfit=fitdata.('param_SCFSpin_Pow');
        eta_vals_fit=fitdata.('eta_vals');
        eta_vals_FS = xyfit_FS.('eta_vals');
        xi_vals=fitdata.('xi_vals');
        i_N = numel(sqrtN_vals);
        r_vals = rbin{i_N,1};
        T_select=[1:22];
        
        L_vals=sqrtN_vals;
        r_min = 3;
        r_max = 55;
        
        
        T_max = 3.5;
        FSplot_min = 1.2;
        FSplot_max = 2.1;
        FS_Tstep = .5;
        T_offset = .03; % For text in inset
        
    end
    
    for i_N = 1:numel(sqrtN_vals)
        for i_T = 1:numel(T_vals)
            chimparq_av{i_N,i_T} = chimparq_av{i_N,i_T} / sqrtN_vals(i_N)^2;
            chimperpq_av{i_N,i_T} = chimperpq_av{i_N,i_T} / sqrtN_vals(i_N)^2;
            chimxq_av{i_N,i_T} = chimxq_av{i_N,i_T} / sqrtN_vals(i_N)^2;
            chimyq_av{i_N,i_T} = chimyq_av{i_N,i_T} / sqrtN_vals(i_N)^2;
            chiteq_av{i_N,i_T} = chiteq_av{i_N,i_T} / sqrtN_vals(i_N)^2;
            chiwq_av{i_N,i_T} = chiwq_av{i_N,i_T} / sqrtN_vals(i_N)^2;
        end
    end


%     [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% Plot qL vs chi
    N_sqrtN = numel(sqrtN_vals);
    c_map = linspecer(numel(sqrtN_vals));
    for i_plot_type = 1:numel(plot_type)
        plot_type_cur = plot_type{i_plot_type};
        for index_T = 1:numel(T_select)
            i_T = T_select(index_T);
            T = T_vals(i_T);
            figure(100 * i_model + 1000 * i_plot_type + i_T);
            for i_N = 1:numel(sqrtN_vals)
                sqrtN = sqrtN_vals(i_N);
                L = L_vals(i_N);
                if (eta_vals_FS(i_T) ~= 0)
                    eta = eta_vals_FS(i_T);
                else
                    eta = eta_vals_fit(end,i_T);
                end
                dispname=sprintf('$N = (%d)^2$', sqrtN);

                if (plot_type_cur == "m")
                    chi_cur = (chimxq_av{i_N,i_T} + chimyq_av{i_N,i_T});
                    hYLabel = ylabel('$\chi_{m}(q)$','interpreter','latex');
                elseif (plot_type_cur == "mpar")
                    chi_cur = chimparq_av{i_N,i_T};
                    hYLabel = ylabel('$\chi_{m\parallel}(q)$','interpreter','latex');
                elseif (plot_type_cur == "mperp")
                    chi_cur = chimperpq_av{i_N,i_T};
                    hYLabel = ylabel('$\chi_{m\perp}(q)$','interpreter','latex');
                elseif (plot_type_cur == "w")
                    chi_cur = chiwq_av{i_N,i_T};
                    hYLabel = ylabel('$\chi_{w}(q)$','interpreter','latex');
                elseif (plot_type_cur == "te")
                    chi_cur = chiteq_av{i_N,i_T};
                    hYLabel = ylabel('$\chi_{\theta}(q)$','interpreter','latex');
                end
                q_vals=qbin{i_N,i_T};
                L_exponent = 0;

                plot(q_vals * L,L^(L_exponent)*chi_cur, ...
                    'Color',c_map(i_N,:),'DisplayName',dispname, ...
                    'LineStyle', ':', ...
                    'LineWidth',2.5,...
                    'Marker', '^', 'MarkerSize', 10);
    %                 'Marker', '^', 'MarkerSize', 4, ...
    %                 'LineWidth',1.5);
                hold on;
            end

            h_axis = gca;
            set(h_axis,'xscale','lin');
    %             xlim([.5,Inf]);
    %             if (plot_type_cur ~= "mparmperp") 
    %                 ylim([-Inf,1.2]);
    %             end
            hXLabel = xlabel('$qL$','interpreter','latex');
    %             hYLabel = ylabel('$C_{m}(q,t)$','interpreter','latex');



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Legend, axes etc
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hLegend = legend('Location', 'Northeast','interpreter','latex','FontSize', 10,...
                'NumColumns',1);

            hTitle = title(curtitle);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Font
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            set(h_axis, 'FontName', 'cmr12')
            set([hXLabel, hYLabel], 'FontName', 'cmr12')
            set(h_axis, 'FontSize', 14)
            set([hXLabel, hYLabel], 'FontSize', 14)
            set(hTitle, 'FontSize', 14, 'FontWeight' , 'bold')


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Adjust axes properties
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
                'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
                'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
                'LineWidth', .5,'Xscale', 'log','Yscale', 'log')

            annotation_str = {sprintf('$T = %.3f$',T),sprintf('$\\eta = %.2f$',eta)};
            dim=[.15 .1 .1 .2];
            annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
                'interpreter','latex',...
                'VerticalAlignment','middle', 'HorizontalAlignment','left',...
                'Color','Red','FontSize', 10,...
                'BackgroundColor','white');
%                 'Units','normalized',...
%                 'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
%                 'interpreter','latex',...
%                 'Color','Red','FontSize', 10);
%             text(0.05,.1,sprintf('$T = %.3f$',T),...
%                 'Units','normalized',...
%                 'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
%                 'interpreter','latex',...
%                 'Color','Red','FontSize', 10);
% 
%             text(0.05,.05,sprintf('$\\eta = %.2f$',eta),...
%                 'Units','normalized',...
%                 'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
%                 'interpreter','latex',...
%                 'Color','Red','FontSize', 10);

            figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/chi_plots/%s_Chi_%s_FS_T_%.3f',curmodel,plot_type_cur,T);
            fprintf('Creating figure %s\n',figname)
            if(saveswitch == 1)
               exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
               exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
            end
        end
    end
    
    
    
    
    
    %% Plot qL vs chi L^{-(2-\eta)}
    N_sqrtN = numel(sqrtN_vals);
    c_map = linspecer(numel(sqrtN_vals));
    for i_plot_type = 1:numel(plot_type)
        plot_type_cur = plot_type{i_plot_type};
        for index_T = 1:numel(T_select)
            i_T = T_select(index_T);
            T = T_vals(i_T);
            figure(100 * i_model + 1000 * i_plot_type + i_T + 1e5);
            for i_N = 1:numel(sqrtN_vals)
                sqrtN = sqrtN_vals(i_N);
                L = L_vals(i_N);
                if (eta_vals_FS(i_T) ~= 0)
                    eta = eta_vals_FS(i_T);
                else
                    eta = eta_vals_fit(end,i_T);
                end
                dispname=sprintf('$N = (%d)^2$', sqrtN);

                if (plot_type_cur == "m")
                    chi_cur = (chimxq_av{i_N,i_T} + chimyq_av{i_N,i_T})/2;
                    hYLabel = ylabel('$L^{-(2-\eta)}\chi_{m}(q)$','interpreter','latex');
                elseif (plot_type_cur == "mpar")
                    chi_cur = chimparq_av{i_N,i_T};
                    hYLabel = ylabel('$L^{-(2-\eta)}\chi_{m\parallel}(q)$','interpreter','latex');
                elseif (plot_type_cur == "mperp")
                    chi_cur = chimperpq_av{i_N,i_T};
                    hYLabel = ylabel('$L^{-(2-\eta)}\chi_{m\perp}(q)$','interpreter','latex');
                elseif (plot_type_cur == "w")
                    chi_cur = chiwq_av{i_N,i_T};
                    hYLabel = ylabel('$L^{-(2-\eta)}\chi_{w}(q)$','interpreter','latex');
                elseif (plot_type_cur == "te")
                    chi_cur = chiteq_av{i_N,i_T};
                    hYLabel = ylabel('$L^{-(2-\eta)}\chi_{\theta}(q)$','interpreter','latex');
                end
                q_vals=qbin{i_N,i_T};
                L_exponent = -(2-eta);

                plot(q_vals * L,L^(L_exponent)*chi_cur, ...
                    'Color',c_map(i_N,:),'DisplayName',dispname, ...
                    'LineStyle', ':', ...
                    'LineWidth',2.5,...
                    'Marker', '^', 'MarkerSize', 10);
    %                 'Marker', '^', 'MarkerSize', 4, ...
    %                 'LineWidth',1.5);
                hold on;
            end

            h_axis = gca;
            set(h_axis,'xscale','lin');
    %             xlim([.5,Inf]);
    %             if (plot_type_cur ~= "mparmperp") 
    %                 ylim([-Inf,1.2]);
    %             end
            hXLabel = xlabel('$qL$','interpreter','latex');
    %             hYLabel = ylabel('$C_{m}(q,t)$','interpreter','latex');



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Legend, axes etc
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hLegend = legend('Location', 'Northeast','interpreter','latex','FontSize', 10,...
                'NumColumns',1);

            hTitle = title(curtitle);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Font
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            set(h_axis, 'FontName', 'cmr12')
            set([hXLabel, hYLabel], 'FontName', 'cmr12')
            set(h_axis, 'FontSize', 12)
            set([hXLabel, hYLabel], 'FontSize', 12)
            set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Adjust axes properties
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
                'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
                'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
                'LineWidth', .5,'Xscale', 'log','Yscale', 'log')

            text(0.05,.1,sprintf('$T = %.3f$',T),...
                'Units','normalized',...
                'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
                'interpreter','latex',...
                'Color','Red','FontSize', 10);

            text(0.05,.05,sprintf('$\\eta = %.2f$',eta),...
                'Units','normalized',...
                'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
                'interpreter','latex',...
                'Color','Red','FontSize', 10);

            figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/chi_plots/%s_Chi_%s_FS_L2-eta_T_%.3f',curmodel,plot_type_cur,T);
            fprintf('Creating figure %s\n',figname)
            if(saveswitch == 1)
               exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
               exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
            end
        end
    end
end
    