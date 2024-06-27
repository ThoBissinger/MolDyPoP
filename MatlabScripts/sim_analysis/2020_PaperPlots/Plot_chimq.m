clear all
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;

mxydata=load('mxy/rho_3.00_dynamics_LinearTime.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat'); mxyfit_FS=load('mxy/rho_3.00_DataCollapse_SCF.mat');
xydata=load('xy/lf0_eq.mat'); xyfit=load('xy/lf0_CritExpFit.mat');


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
        
        plot_type={"mpar", "mperp", "m", "w", "te"};
        
        expfit=fitdata.('param_SCFSpin_Exp');
        powfit=fitdata.('param_SCFSpin_Pow');
        eta_vals=fitdata.('eta_vals');
        xi_vals=fitdata.('xi_vals');
        
        i_N = numel(sqrtN_vals);
        r_vals = rbin{i_N,1};
        T_select=[9:2:13,14:21];
%         T_select=1:numel(T_vals);
%         T_select=[2,3,4,5,7];
        
        L_vals=[9.25,18.5,37,74,148];
        r_min = 2;
        r_max = 60;
        
        T_max = .4;
        FSplot_min = .13;
        FSplot_max = .23;
        FS_Tstep = .05;
        T_offset = .0025; % For text in inset
%         labels=['N = 2^{10}', 'N = 2^{12}', 'N = 2^{14}', 'N = 2^{16}'];
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
        
        expfit=fitdata.('param_SCFSpin_Exp');
        powfit=fitdata.('param_SCFSpin_Pow');
        eta_vals=fitdata.('eta_vals');
        xi_vals=fitdata.('xi_vals');
        
        i_N = numel(sqrtN_vals);
        r_vals = rbin{i_N,1};
        T_select=[17:3:38];
        
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
    %% Create basic plot
     %% Plot
    N_sqrtN = numel(sqrtN_vals);
    c_map = linspecer(numel(T_select));
    for i_plot_type = 1:numel(plot_type)
        plot_type_cur = plot_type{i_plot_type};
        figure(50 * i_model + 1000 * i_plot_type)
        hold on;
        for index_T = 1:numel(T_select)
            i_T = T_select(index_T);
            T = T_vals(i_T);
            
            dispname=sprintf('T = %.3f', T);

            if (plot_type_cur == "m")
                chi_cur = (chimxq_av{i_N,i_T} + chimyq_av{i_N,i_T})/2;
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
%                 TCF_Spin_cur = chi_cur(i_q:length(q_vals):end);

%                 if (plot_type_cur ~= "mparmperp") % No division by value at 0 for this one, because it should be zero anyway.
% %                 if (orientation_cur ~= "mparmperp" && orientation_cur ~= "mparmpar") 
%                     TCF_Spin_cur = TCF_Spin_cur / TCF_Spin_cur(1);
%                 end
%             gxx_cur = gxx{i_N,i_T};
%                 t_vals=TCF_times{i_N,i_T};
            q_vals=qbin{i_N,i_T};
%             TCF_Spin_cur = gxx_cur(i_q:length(q_vals):end);
%                 fitfunc_DO=@(c) TCF_Spin_cur(1)*exp(-c(1) * t_vals/2) .* cos(c(2) * t_vals - c(3))/cos(c(3));
    %         c_forq = coeffs_cur(3*(i_q - 1) + (1:3));

            plot(q_vals,chi_cur, ...
                'Color',c_map(index_T,:),'DisplayName',dispname, ...
                'LineStyle', '-', ...
                'LineWidth',1.8);
%                 'Marker', '^', 'MarkerSize', 4, ...
%                 'LineWidth',1.5);
            hold on;
    %         legend show
    %         legend('Location','best')
        end

        h_axis = gca;
        set(h_axis,'xscale','lin');
%             xlim([.5,Inf]);
%             if (plot_type_cur ~= "mparmperp") 
%                 ylim([-Inf,1.2]);
%             end
        hXLabel = xlabel('$q$','interpreter','latex');
%             hYLabel = ylabel('$C_{m}(q,t)$','interpreter','latex');



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Legend, axes etc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hLegend = legend('Location', 'Best','interpreter','latex','FontSize', 10,...
            'NumColumns',3);

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
            'LineWidth', .5,'Xscale', 'lin','Yscale', 'log')

%             text(0.05,.45,sprintf('$N = %d^2$',sqrtN_vals(i_N)),...
%                 'Units','normalized',...
%                 'VerticalAlignment','top', 'HorizontalAlignment','left',...
%                 'interpreter','latex',...
%                 'Color','Red','FontSize', 10);

%             text(0.05,.40,sprintf('$q = %.3f$',q),...
%                 'Units','normalized',...
%                 'VerticalAlignment','top', 'HorizontalAlignment','left',...
%                 'interpreter','latex',...
%                 'Color','Red','FontSize', 10);

        figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/%s_Chi_%s_sqrtN_%d',curmodel,plot_type_cur,sqrtN_vals(i_N));
        fprintf('Creating figure %s\n',figname)
        if(saveswitch == 1)
           exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
    end
    
   
end
    