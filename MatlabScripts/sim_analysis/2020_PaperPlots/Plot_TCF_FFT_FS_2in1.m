clear
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/TCF/FFT_FS_2in1',fig_base);

section_select = 1;

mxydata=matfile(mxydata_AdjustedTime_name);
fmxydata=matfile(mxydata_AdjustedTime_name);
        
i_q=9;
n_q = 3;

T_select=[1:7];
sqrtN_select=[2:4];
xlim_max_factor = 2;
xlim_max = 20;

L_vals=[9.25,18.5,37,74,148];

res_factor=3;
res_function=@(t,tau) resolution_Laplace_pleateau(t,tau,tau);

i_N_start = 2;
i_N_max = 4;
n_trunc = 4; % For NF Psi
plot_type=["mperp", "m"];
% All the data for mxy model
if (ismember("m",plot_type))
    mxy_gxx=mxydata.('gxx');
    mxy_gyy=mxydata.('gyy');
    mxy_chimxq_av=mxydata.('chimxq_av');
    mxy_chimyq_av=mxydata.('chimyq_av');
end
if (ismember("mpar",plot_type))
    mxy_gmparmpar=mxydata.('gmparmpar');
    mxy_chimparq_av=mxydata.chimparq_av;
end
if (ismember("mperp",plot_type))
    mxy_gmperpmperp=mxydata.('gmperpmperp');
    mxy_chimperpq_av=mxydata.chimperpq_av;
end
if (ismember("t",plot_type))
    mxy_gtt=mxydata.('gtt');
    % TODO chiteq. Missing in current dataset
end
if (ismember("w",plot_type))
    mxy_gww=mxydata.('gww');
    mxy_chiwq_av=mxydata.('chiwq_av');
end
    
mxy_qbin=mxydata.('qbin');
mxy_averaging_times=mxydata.('averaging_times');
mxy_absM_av=mxydata.('absM_av');

mxy_sqrtN_vals=mxydata.('sqrtN_vals');
mxy_T_vals=mxydata.('T_vals');

% All the data for fmxy model
if (ismember("m",plot_type))
    fmxy_gxx=fmxydata.('gxx');
    fmxy_gyy=fmxydata.('gyy');
    fmxy_chimxq_av=fmxydata.('chimxq_av');
    fmxy_chimyq_av=fmxydata.('chimyq_av');
end
if (ismember("mpar",plot_type))
    fmxy_gmparmpar=fmxydata.('gmparmpar');
    fmxy_chimparq_av=fmxydata.chimparq_av;
end
if (ismember("mperp",plot_type))
    fmxy_gmperpmperp=fmxydata.('gmperpmperp');
    fmxy_chimperpq_av=fmxydata.chimperpq_av;
end
if (ismember("t",plot_type))
    fmxy_gtt=fmxydata.('gtt');
    % TODO chiteq. Missing in current dataset
end
if (ismember("w",plot_type))
    fmxy_gww=fmxydata.('gww');
    fmxy_chiwq_av=fmxydata.('chiwq_av');
end
    
fmxy_qbin=fmxydata.('qbin');
fmxy_averaging_times=fmxydata.('averaging_times');
fmxy_absM_av=fmxydata.('absM_av');

fmxy_sqrtN_vals=fmxydata.('sqrtN_vals');
fmxy_T_vals=fmxydata.('T_vals');

z=1;
for i_plot_type = 1:numel(plot_type)
    figure
    plot_type_cur = plot_type(i_plot_type);
    %% Section 1 MXY Model
    subplot(2,1,2);
    N_N=numel(sqrtN_select);
    N_T = numel(T_select);
    c_map = linspecer(N_N);
    max_ft_vec=zeros(1,N_N);
    S_q_vec=zeros(1,N_N);
    om_max_vec=zeros(1,N_N);
    for ind_T = 1:N_T
        i_T = T_select(ind_T);
        T = mxy_T_vals(i_T);
        figure
        hold on;
        for ind_N = 1:N_N
            i_N = sqrtN_select(ind_N);
%             i_T = pairs(i_pair,2);
%             i_q = pairs(i_pair,3);
%             n_q = i_q;
            L = L_vals(i_N);
    
            sqrtN = mxy_sqrtN_vals(i_N);
            
            q_vals = mxy_qbin{i_N,i_T};
            q=q_vals(i_q);

            t=mxy_averaging_times{i_N,i_T};
            N_t = numel(t);
            N_q = numel(q_vals);
            q_indices=i_q:N_q:N_q*N_t;    
            if (plot_type_cur == "m")
                cf = mxy_gxx{i_N,i_T}(q_indices) + mxy_gyy{i_N,i_T}(q_indices);
            elseif (plot_type_cur == "mpar")
                cf = mxy_gmparmpar{i_N,i_T}(q_indices);
            elseif (plot_type_cur == "mperp")
                cf = mxy_gmperpmperp{i_N,i_T}(q_indices);
            elseif (plot_type_cur == "te")
                cf = mxy_gtt{i_N,i_T}(q_indices);
            elseif (plot_type_cur == "w")
                cf = mxy_gww{i_N,i_T}(q_indices);
            end
            c = fit_DampedOscillator_RealSpace(t,cf,n_period,weightexp,'omega_1');
            gamma_cur = c(1);
            omega_1_cur = c(2);
            omega_0_cur = sqrt(c(2)^2 - gamma_cur^2/4);
    
            tau=res_factor/gamma_cur;
            resolution_vals=res_function(t,tau);
                
            [ft_vals,om_vals]=FT_correlation(t, cf .* resolution_vals, 1e5);
            ft_vals=real(ft_vals);
            S_q = cf(1);
            [ft_max,i_ft_max] = max(ft_vals);
            
            max_ft_vec(ind_N) = ft_max;
            S_q_vec(ind_N) = S_q;
            om_max_vec(ind_N) = abs(om_vals(i_ft_max));
    
            dispname = sprintf('$N = (%d)^2$',sqrtN);
%             loglog(om_vals / om_vals(i_ft_max), real(ft_vals)/ L^z / S_q,...
            loglog(om_vals * L^z, real(ft_vals)/ L^z / S_q,...
                'Color',c_map(ind_N,:),'DisplayName',dispname, ...
                'LineStyle', '-', 'LineWidth',2, ...
                'Marker', 'none', 'MarkerSize', 3);
            hold on;
        end
        absM_vec=cell2mat(absM_av(:,i_T));
        eta_fitob=fit_eta_Magnetization_FS(absM_vec,L_vals);
        eta=eta_fitob.eta;
%         y_vals=logspace(log10(xlim_vec(1)),log10(xlim_vec(2)),1e2);
%         NF_Psi = NelsonFisher_Psi(y_vals,eta,n_trunc);
% %                 loglog(om_vals,2*pi^2*eta^2*om_vals.^(-3-eta),...
%         loglog(y_vals,NF_Psi,...
%                 'Color','Black','DisplayName','NF', ...
%                 'LineStyle', '--', 'LineWidth',2, ...
%                 'Marker', 'none', 'MarkerSize', 5);

        if xlim_max_factor * om_max_vec(end)*L(end)^z > xlim_max
            xlim([0 xlim_max_factor*om_max_vec(end)*L(end)^z]);
        else
            xlim([0 xlim_max]);
        end
        ylim([0 Inf]);
%         ylim([.3*min(NF_Psi),10*max(NF_Psi)]);
        h_axis = gca;
        hXLabel = xlabel('$\omega L^z$','interpreter','latex');
        if (plot_type_cur == "m")
            hYLabel = ylabel('$S_{m}(q,\omega) / \chi_{m}(q) L^z$','interpreter','latex');
        elseif (plot_type_cur == "mpar")
            hYLabel = ylabel('$S_{m\parallel}(q,\omega) / \chi_{m\parallel}(q) L^z$','interpreter','latex');
        elseif (plot_type_cur == "mperp")
            hYLabel = ylabel('$S_{m\perp}(q,\omega) / \chi_{m\perp}(q) L^z$','interpreter','latex');
        elseif (plot_type_cur == "te")
            hYLabel = ylabel('$S_{\theta}(q,\omega) / \chi_{\theta}(q) L^z$','interpreter','latex');
        elseif (plot_type_cur == "w")
            hYLabel = ylabel('$S_{w}(q,\omega) / \chi_{w}(q) L^z$','interpreter','latex');
        end
        
        



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Legend, axes etc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hLegend = legend('Location', 'Northeast','interpreter','latex',...
            'NumColumns',1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Font
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set_fonts_default

        set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
            'LineWidth', .5)

        annotation_str = {model_str,sprintf('$T = %.3f$',T),sprintf('$q = %d \\cdot 2\\pi / L$',n_q),... sprintf('$n_q = %d$',i_q)
            sprintf('$\\eta = %.2f$',eta),sprintf('$z = %.2f$',z)};
        dim=[.68 .45 .1 .2];
        annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
            'interpreter','latex',...
            'VerticalAlignment','middle', 'HorizontalAlignment','left',...
            'Color','black','FontSize', fontsize_annotation,...
            'BackgroundColor','white');

        set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

%                 figname=sprintf('%s/HighFreq/%s/%s_S_FS_n_q_%d_z_%.2f_T_%.3f_Nstart_%d_resfac_%d',basedir,plot_type_cur,curmodel,n_q,z,T,i_N_start,res_factor);
        figname=sprintf('%s/%s_%s_TCF_FFT_FS_T_%.3f_nq_%d',basedir,plot_type_cur,curmodel,T,n_q);
%                 figname=sprintf('%s/SbyL_vsomL/%s_S_%s_FS_n_q_%d_z_%.2f_tau_%de3_T_%.3f',basedir,curmodel,plot_type_cur,n_q,z,tau_laplace/1e3,T);
        fprintf('Creating figure %s\n',figname)
        if(saveswitch == 1)
           exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');

           set(gca,'YScale','log')
           figname=sprintf('%s/%s_%s_TCF_FFT_FS_logy_T_%.3f_nq_%d',basedir,plot_type_cur,curmodel,T,n_q);
           ylim([1e-2 Inf])
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end

    end        
end
