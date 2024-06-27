clear all
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=0;

% mxydata=load('mxy/rho_3.00_integ.mat');
% mxydata=load('mxy/rho_3.00_dynamics_better_q.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat');
mxydata=load('mxy/rho_3.00_dynamics_LinearTime.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat');mxyfit_FS=load('mxy/rho_3.00_DataCollapse_SCF.mat');
xydata=load('xy/lf0_qreduced.mat'); xyfit=load('xy/lf0_CritExpFit.mat');


for i_model = 1:1

    % figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=mxydata;
        fitdata=mxyfit_FS;
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        gxx=data.('gxx');
        gyy=data.('gyy');
        gmparmpar=data.('gmparmpar');
        gmperpmperp=data.('gmperpmperp');
        gmparmperp=data.('gmparmperp');
        TCF_times=data.('TCF_times');
        qbin = data.('qbin');
        ACF_q0_M = data.('ACF_q0_M');
                
        i_N = numel(sqrtN_vals);
        q_vals = qbin{i_N,1};
        T_select=[9:2:13,14:21];
        T_select=1:numel(T_vals);
%         T_select=[2,3,4,5,7];
        
        orientations={"mm"};
%         orientations={"mm","mparmpar", "mperpmperp","mparmperp"};
        
        L_vals=[9.25,18.5,37,74,148];
        r_min = 6;
        r_max = 35;
    
        q_select=0:8;
%         q_select=0:10;
        
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
        r_min = 8;
        r_max = 55;
        
        
        T_max = 3.5;
        FSplot_min = 1.2;
        FSplot_max = 2.1;
        FS_Tstep = .5;
        T_offset = .03; % For text in inset
        
    end


%     [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% (gmxx + gmyy) / 2
    N_sqrtN = numel(sqrtN_vals);
    c_map = linspecer(numel(T_select));
    for i_orientation = 1:numel(orientations);
        orientation_cur = orientations{i_orientation};
        for i_q = q_select
    %         i_q = q_select(index_q);
            if (i_q == 0)
                q = 0;
                if (orientation_cur ~= "mm")
                    continue
                end
            else
                q = q_vals(i_q);
            end
            figure(i_q+1 + 50 * i_model + 1000 * i_orientation)
        %     coeffs = cell2mat(param_TCFSpin_omega_1_q_DO(i_N,:)');
            for index_T = 1:length(T_select)
        %         fig1=figure(10*i_model + i);
        %         hold on
                i_T = T_select(index_T);

                T = T_vals(i_T);
                dispname=sprintf('T = %.3f', T);
                if ( q == 0 )
                    TCF_Spin_cur = ACF_q0_M{i_N,i_T};
                    hYLabel = ylabel('$C_{m}(q,t)$','interpreter','latex');
    %                 TCF_Spin_cur = gxx_cur(i_q:length(q_vals):end);
                else
                    if (orientation_cur == "mm")
                        g_cur = (gxx{i_N,i_T} + gyy{i_N,i_T})/2;
                        hYLabel = ylabel('$C_{m}(q,t)$','interpreter','latex');
                    elseif (orientation_cur == "mparmpar")
                        g_cur = gmparmpar{i_N,i_T};
                        hYLabel = ylabel('$C_{m\parallel}(q,t)$','interpreter','latex');
                    elseif (orientation_cur == "mperpmperp")
                        g_cur = gmperpmperp{i_N,i_T};
                        hYLabel = ylabel('$C_{m\perp}(q,t)$','interpreter','latex');
                    elseif (orientation_cur == "mparmperp")
                        g_cur = gmparmperp{i_N,i_T};
                        hYLabel = ylabel('$C_{m\parallel\perp}(q,t)$','interpreter','latex');
                    end
                    TCF_Spin_cur = g_cur(i_q:length(q_vals):end);
                end
                
                if (orientation_cur ~= "mparmperp") % No division by value at 0 for this one, because it should be zero anyway.
%                 if (orientation_cur ~= "mparmperp" && orientation_cur ~= "mparmpar") 
                    TCF_Spin_cur = TCF_Spin_cur / TCF_Spin_cur(1);
                end
    %             gxx_cur = gxx{i_N,i_T};
                t_vals=TCF_times{i_N,i_T};
    %             TCF_Spin_cur = gxx_cur(i_q:length(q_vals):end);
                fitfunc_DO=@(c) TCF_Spin_cur(1)*exp(-c(1) * t_vals/2) .* cos(c(2) * t_vals - c(3))/cos(c(3));
        %         c_forq = coeffs_cur(3*(i_q - 1) + (1:3));

                plot(t_vals,real(TCF_Spin_cur ), ...
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
            set(h_axis,'xscale','log');
            xlim([.5,Inf]);
            if (orientation_cur ~= "mparmperp") 
                ylim([-Inf,1.2]);
            end
            hXLabel = xlabel('$t$','interpreter','latex');
%             hYLabel = ylabel('$C_{m}(q,t)$','interpreter','latex');



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Legend, axes etc
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hLegend = legend('Location', 'SouthWest','interpreter','latex','FontSize', 10,...
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
                'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')

            text(0.05,.45,sprintf('$N = %d^2$',sqrtN_vals(i_N)),...
                'Units','normalized',...
                'VerticalAlignment','top', 'HorizontalAlignment','left',...
                'interpreter','latex',...
                'Color','Red','FontSize', 10);

            text(0.05,.40,sprintf('$q = %.3f$',q),...
                'Units','normalized',...
                'VerticalAlignment','top', 'HorizontalAlignment','left',...
                'interpreter','latex',...
                'Color','Red','FontSize', 10);

            figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/%s_SpinTCF_%s_sqrtN_%d_q_%.3f',curmodel,orientation_cur,sqrtN_vals(i_N),q);
            fprintf('Creating figure %s\n',figname)
            if(saveswitch == 1)
               exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
               exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
            end
        end
    end
    
    
end
    