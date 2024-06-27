clear all
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;

% mxydata=load('mxy/rho_3.00_integ.mat');
mxydata=load('mxy/rho_3.00_dynamics_LinearTime.mat');
mxyfit=load('mxy/rho_3.00_TCF_q_fit.mat');

xydata=load('xy/xy_dynamics_LinearTime.mat');
xyfit=load('xy/lf0_CritExpFit.mat');


for i_model = 1:2

    % figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=mxydata;
        fitdata=mxyfit;
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        
        gtt=data.('gtt');
        TCF_times=data.('TCF_times');
        qbin = data.('qbin');
        
        i_N = numel(sqrtN_vals);
        q_vals = qbin{i_N,1};
        T_select=[9:2:13,14:21];
        T_select=1:numel(T_vals);
        T_select=1:20;
        
        L_vals=[9.25,18.5,37,74,148];
        r_min = 6;
        r_max = 35;
    
        q_select=1:6;
        

    else
        curmodel="xy";
        curtitle="SXY model";
        
        data=xydata;
        fitdata=xyfit;

        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        
        gtt=data.('gtt');
        TCF_times=data.('TCF_times');
        qbin = data.('qbin');
        
        i_N = numel(sqrtN_vals);
        q_vals = qbin{i_N,1};
%         T_select=[9:2:13,14:21];
%         T_select=1:numel(T_vals);
        T_select=1:20;
        
        L_vals=[9.25,18.5,37,74,148];
        r_min = 6;
        r_max = 35;
    
        q_select=1:6;
        
        
    end


%     [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% Create basic plot
    N_sqrtN = numel(sqrtN_vals);
    c_map = linspecer(numel(q_select));
    
    for img_index = 1:length(T_select)
        figure(img_index)
        i_T = T_select(img_index);
        T = T_vals(i_T);
        
        gtt_cur = gtt{i_N,i_T};
        
        ymax=0;
        ymin=0;
        for i_q = 1:numel(q_select)
            q = q_vals(i_q);
            dispname=sprintf('q = %.3f', q);
    %     coeffs = cell2mat(param_TCFSpin_omega_1_q_DO(i_N,:)');
        

            
            t_vals=TCF_times{i_N,i_T};
            TCF_cur = gtt_cur(i_q:length(q_vals):end);

            plot(q * t_vals,q^2 * real(TCF_cur),'.',...
                'Color',c_map(i_q,:),'DisplayName',dispname,...
                'LineStyle', '-', ...
                'LineWidth',1.5);
%                 'Marker', '^', 'MarkerSize', 4, ...
%                 'LineWidth',1.5);
            hold on;
            
            ymax=max(ymax,q^2 * max(abs(TCF_cur)));
            ymin=min(ymin,q^2 * min(real(TCF_cur)));
    %         legend show
    %         legend('Location','best')
        end

        h_axis = gca;
        set(h_axis,'xscale','log');

        hXLabel = xlabel('$qt$','interpreter','latex');
        hYLabel = ylabel('$q^2 C_{\theta\theta}(q,t)$','interpreter','latex');



%         ylim([min(1.1*ymin,-.4*ymax),1.3*ymax]);
        xlim([0,500]);

        %% Legend, axes etc
        hLegend = legend('Location', 'Best','interpreter','latex','FontSize', 10,...
            'NumColumns',2);

        hTitle = title(curtitle);

        %% Font
        set(h_axis, 'FontName', 'cmr12')
        set([hXLabel, hYLabel], 'FontName', 'cmr12')
        set(h_axis, 'FontSize', 12)
        set([hXLabel, hYLabel], 'FontSize', 12)
        set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


        %% Adjust axes properties
        set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
            'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')

        text(0.95,.95,sprintf('$N = %d^2$',sqrtN_vals(i_N)),...
            'Units','normalized',...
            'VerticalAlignment','top', 'HorizontalAlignment','right',...
            'interpreter','latex',...
            'Color','Blue','FontSize', 10);

        text(0.95,.90,sprintf('$T = %.3f$',T),...
            'Units','normalized',...
            'VerticalAlignment','top', 'HorizontalAlignment','right',...
            'interpreter','latex',...
            'Color','Blue','FontSize', 10);

        figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/TCF/tete/%s_teteTCF_sqrtN_%d_T_%.3f',curmodel,sqrtN_vals(i_N),T);
        fprintf('Creating figure %s\n',figname)
        if(saveswitch == 1)
           exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end

    end
    
    
%     for img_index = 1:length(T_select)
%         figure(img_index+100)
%         i_T = T_select(img_index);
%         T = T_vals(i_T);
%         
%         gtt_cur = gxx{i_N,i_T};
%         
%         curcoeff = param_TCFSpin_omega_1_q_DO{i_N,i_T};
%         
%         ymax=0;
%         ymin=0;
%         for i_q = 1:numel(q_select)
%             q = q_vals(i_q);
%             dispname=sprintf('q = %.3f', q);
%         coeffs = cell2mat(param_TCFSpin_omega_1_q_DO(i_N,:)');
%         
% 
%             
%             t_vals=TCF_times{i_N,i_T};
%             TCF_cur = gtt_cur(i_q:length(q_vals):end);
% 
%             plot(t_vals(1:4:end),real(TCF_cur(1:4:end)),...
%                 'Color',c_map(i_q,:),...
%                 'LineStyle', 'none', 'HandleVisibility', 'off',... % 'IconDisplayStyle', 'off', ...
%                 'Marker', 'o', 'MarkerSize', 5, 'LineWidth',1.8);
%             hold on;
%                         
%             plot(t_vals,fitfunc_DO(t_vals,real(TCF_cur(1)),curcoeff(3*(i_q-1)+(1:3))),...
%                 'Color',c_map(i_q,:),'DisplayName',dispname,...
%                 'LineStyle', '-', ...
%                 'LineWidth',1.8);
%             
%             ymax=max(ymax,max(abs(TCF_cur)));
%             ymin=min(ymin,min(real(TCF_cur)));
%             legend show
%             legend('Location','best')
%         end
% 
%         h_axis = gca;
%         set(h_axis,'xscale','log');
% 
%         hXLabel = xlabel('$t$','interpreter','latex');
%         hYLabel = ylabel('$C_{m}(q,t)$','interpreter','latex');
% 
% 
% 
%         ylim([min(1.3*ymin,-.4*ymax),1.3*ymax]);
%         xlim([t_vals(1),t_vals(end)]);
% 
%         % Legend, axes etc
%         hLegend = legend('Location', 'SouthWest','interpreter','latex','FontSize', 10,...
%             'NumColumns',2);
% 
%         hTitle = title(curtitle);
% 
%         % Font
%         set(h_axis, 'FontName', 'cmr12')
%         set([hXLabel, hYLabel], 'FontName', 'cmr12')
%         set(h_axis, 'FontSize', 12)
%         set([hXLabel, hYLabel], 'FontSize', 12)
%         set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
% 
% 
%         % Adjust axes properties
%         set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
%             'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
%             'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
%             'LineWidth', .5,'Xscale', 'log','Yscale', 'lin')
% 
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
% 
%         figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/%s_SpinTCF_sqrtN_%d_T_%.3f',curmodel,sqrtN_vals(i_N),T);
%         fprintf('Creating figure %s\n',figname)
%         if(saveswitch == 1)
%            exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
%            exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
%         end
% 
%     end
end
    