%% 0 Initialization
clear
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/TCF/PeakAndLowFreq',fig_base);

i_model = 4;
if (i_model == 1)
    curmodel="mxy";
    curtitle="MXY model";
    modelname="MXY";

%     sqrtN_vals = [16 32 64 128 256];
    sqrtN_vals = [16 32 64 128];
    dir_longtime="/data/scc/thobi/211201_LongerTime/mxy_3.00";
    dir_shorttime="/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00";

    T_vals = [.11 .14 .17 .185];
    T_dirs = {"T_.11" "T_.14" "T_.17" "T_.185"};
%     T_vals = [.17];
%     T_dirs = {"T_.17"};
%     T_vals = [.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52];
%     T_dirs = {"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
    
    sampfilename="samp_Dynamics";
    
    L_vals=[9.25,18.5,37,74,148];
    q_select=[1 2 3 4 5 6 7 8 10 11 13 14 15 18 21 22]; % Ignores multiple counts due to binning
    q_indices=[1 3 9 23];

    i_N_plot =4;
    i_q_plot=9;
    q_select = [1,3,9];
%     i_N_start = 2;
%     i_N_max = 4;

    runmax_vals=[1e3 1e3 500 500];
        
elseif (i_model == 2)
        curmodel="xy";
        curtitle="XY model";

elseif (i_model == 3)
    curmodel="fmxy";
    curtitle="FMXY model";
    modelname="FMXY";

%     sqrtN_vals = [16 32 64 128 256];
    sqrtN_vals = [16 32 64 128];
    dir_longtime="/data/scc/thobi/211201_LongerTime/fmxy";
    dir_shorttime="/data/scc/thobi/220201_ReducedSmapleStepDeltat/fmxy";

    T_vals = [.11 .14 .17 .185];
    T_dirs = {"T_.11" "T_.14" "T_.17" "T_.185"};
    
    sampfilename="samp_Dynamics";
    
    L_vals=[9.25,18.5,37,74,148];
    
    i_N_plot = 2;
    i_q_plot=1;

    runmax_vals=[3000 1500 500 125];
    q_select = [1,3,9];
    
elseif (i_model == 4)
    curmodel="xy_s";
    curtitle="XY S model";
    modelname="XY";

    sqrtN_vals = [16 32 64 128];
    dir_longtime="/data/scc/thobi/211201_LongerTime/xy_s";
    dir_shorttime="/data/scc/thobi/220201_ReducedSmapleStepDeltat/xy_s";

    T_vals = [.85 .91 .95 1.00];
    T_dirs = {"T_.85" "T_.91" "T_.95" "T_1.00"};
    
    sampfilename="samp_Dynamics";
    
    L_vals=sqrtN_vals;
    q_indices=[1 3];
    i_N_start = 3;
    i_N_max = 4;

    i_N_plot = 4;
    i_q_plot=1;

    runmax_vals=[500 500 500 500];

    q_select = [1,3];
    
end
res_factors=3/4*log(runmax_vals);
res_function = resolution_Gauss;
% res_function=@(t,tau) resolution_Laplace_pleateau(t,tau,tau);

% varname='gmperpmperp'; varshort='mperp';
% varname='gxx'; varshort='x';
% varname='gww'; varshort='w';
varname='mperp'; ylabel_name='$S_{m\perp}(q,\omega)/\chi_{m\perp}(q)$';
varname='m'; ylabel_name='$S_{m}(q,\omega)/\chi_{m}(q)$';
varname='w'; ylabel_name='$S_{w}(q,\omega)/\chi_{w}(q)$';


q_index_start = 6;

n_trunc = 4; % For NF Psi
z=1;

N_T = numel(T_vals);
N_N = numel(sqrtN_vals);
N_q = numel(q_select);
    

%% 1 Data collection
% n_N = i_N_max - i_N_start + 1;
q_vals = zeros(N_N,N_T,N_q);
gamma_vals = zeros(N_N,N_T,N_q);
omega_1_vals = zeros(N_N,N_T,N_q);
om_max_vals = zeros(N_N,N_T,N_q);
ft_max_vals = zeros(N_N,N_T,N_q);
ft_zero_vals = zeros(N_N,N_T,N_q);
ft_sidepeak_max = zeros(N_N,N_T,N_q);
om_vals = cell(N_N,N_T,N_q);
ft_mean = cell(N_N,N_T,N_q);
ft_std = cell(N_N,N_T,N_q);
ft_var = cell(N_N,N_T,N_q);
n_period=4;
weightexp=1;
smooth_depth=5*[1,1,1,1];

% ylim_vec=[1e-3 Inf];
% eta=.24;
n_trunc = 5;

spline_factor=0;

xlim_vec=[1e-1 10*max((spline_factor)^(1/3),1)];

c_map=linspecer(N_T);
eta_vals = zeros(1,N_T);
for i_T = 1:N_T
    T_dir=T_dirs{i_T};
    absM_vec=zeros(1,N_N);
    for i_N = 1:N_N
        sqrtN = sqrtN_vals(i_N);
        curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_longtime,sqrtN,T_dir,sampfilename);
        load(curfile,"absM_av");
        absM_vec(i_N) = absM_av;
    end
    eta_fitob = fit_eta_Magnetization_FS(absM_vec(1:N_N),L_vals(1:N_N));
    eta_vals(i_T) = eta_fitob.eta;
end

% i_N = i_N_plot;
% i_q = i_q_plot;
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATING FFT WITH ERRORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    L = L_vals(i_N);

    for i_T = 1:N_T
        T = T_vals(i_T);
        curdir=sprintf('%s/sqrtN_%d/%s',dir_longtime,sqrtN,T_dirs{i_T});
        curfile_collect=sprintf('%s/%s_collect.mat',curdir,sampfilename);
        curfile_av=sprintf('%s/%s.mat',curdir,sampfilename);
        
        load(curfile_av,'averaging_times','qbin');

        for ind_q = 1:N_q
            i_q = q_select(ind_q);
            q_vals(i_N,i_T,ind_q) = qbin(i_q);
            t=averaging_times;
            if varname == "mperp"
                load(curfile_collect,'gmperpmperp_collect');
                cf=real(gmperpmperp_collect(:,i_q:numel(qbin):end));
            elseif varname == "m"
                load(curfile_collect,'gxx_collect','gyy_collect');
                cf=real(gxx_collect(:,i_q:numel(qbin):end)+gyy_collect(:,i_q:numel(qbin):end));
            elseif varname == "w"
                load(curfile_collect,'gww_collect');
                cf=real(gww_collect(:,i_q:numel(qbin):end));
            end
            cf_mean=mean(cf);
            
            c = fit_DampedOscillator_RealSpace(t,cf_mean,n_period,weightexp,'omega_1');
            gamma_vals(i_N,i_T,ind_q) = c(1);
            omega_1_vals(i_N,i_T,ind_q) = c(2);
        
        
            tau=res_factors(i_N)/gamma_vals(i_N,i_T,ind_q);
            res_vals=res_function(t,tau);
        
            [ft_vals,om_vals_cur]=FT_correlation(t, (cf .* res_vals)' /mean(cf(:,1)), 0);
            ft_vals=real(ft_vals)';
            ft_mean{i_N,i_T,ind_q}=mean(real(ft_vals));
            ft_std{i_N,i_T,ind_q}=std(real(ft_vals));
            ft_var{i_N,i_T,ind_q}=var(real(ft_vals));
            om_vals{i_N,i_T,ind_q}=om_vals_cur;
        
            [ft_max,i_max]=max(ft_mean{i_N,i_T,ind_q});
            om_max_vals(i_N,i_T,ind_q) = abs(om_vals_cur(i_max));
            ft_max_vals(i_N,i_T,ind_q) = ft_max;
        
            i_zero = find(om_vals_cur >= 0, 1);
            ft_zero_vals(i_N,i_T,ind_q) = ft_mean{i_N,i_T,ind_q}(i_zero);
        
            i_sidepeak = find(abs(om_vals_cur) < (om_max_vals(i_N,i_T,ind_q) + gamma_vals(i_N,i_T,ind_q))/2);
            ft_sidepeak_max(i_N,i_T,ind_q) = max(ft_mean{i_N,i_T,ind_q}(i_sidepeak));
        end
    end
end

%% 2 The plot
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    L = L_vals(i_N);

    for ind_q = 1:N_q
        i_q = q_select(ind_q);
        figure
        for i_T = 1:N_T
            T = T_vals(i_T);
            q_cur = q_vals(i_N,i_T,ind_q);

            om_vals_cur=om_vals{i_N,i_T,ind_q};
            ft_mean_cur=ft_mean{i_N,i_T,ind_q};
            ft_std_cur=ft_std{i_N,i_T,ind_q};
        
            errorbar(om_vals_cur,ft_mean_cur,ft_std_cur/sqrt(runmax_vals(i_N)),...
                'HandleVisibility','off',...
                'LineWidth',1,...
                'Color',c_map(i_T,:));
            hold on; 
        end
        for i_T = 1:N_T
            T = T_vals(i_T);
            if i_model == 4
                dispname=sprintf('$T = %.2f$',T);
            else
                dispname=sprintf('$T = %.3f$',T);
            end
        
            om_vals_cur=om_vals{i_N,i_T,ind_q};
            ft_mean_cur=ft_mean{i_N,i_T,ind_q};
            ft_std_cur=ft_std{i_N,i_T,ind_q};
            plot(om_vals_cur,smoothen(ft_mean_cur,1,1),...
                '-',...
                'HandleVisibility','off',...
                'LineWidth',1,...
                'Color','black');
            plot(NaN,NaN,...
                '-',...
                'DisplayName',dispname,...
                'LineWidth',2,...
                'Color',c_map(i_T,:));        
        end


        c_vals = omega_1_vals./q_vals;
        
        eta_cur = eta_vals(i_T);
        
        xlim([0 3*max(om_max_vals(i_N,:,ind_q))]);
        ylim([0 Inf]);
        hLegend=legend('Location','NorthEast','Interpreter','latex',...
            'NumColumns',1);
        set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
            'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')
        
        hXLabel = xlabel('$\omega$','interpreter','latex');
        hYLabel = ylabel(ylabel_name,'interpreter','latex');
        h_axis = gca;
        
        % annotation_str = {modelname,sprintf('$N = (%d)^2$',sqrtN)};
        % annotation_str = {modelname,sprintf('$n_{\\textrm{spline}} = %d$',spline_factor)};
        % annotation_str = {modelname,sprintf('$T = %.3g$',T),sprintf('$q = %.3g$',q_vals(end))};
        annotation_str = {modelname,sprintf('$N = (%d)^2$',sqrtN),sprintf('$q = %.3f$',q_cur)};
        
        % dim=[.16 .45 .1 .2]; % Below left legend
%         dim=[.18 .4 .1 .2]; % Below left legend
        dim=[.595 .4 .1 .2]; % Below of right legend
        % dim=[.7 .7 .1 .2]; % top right corner
        h_annotation = annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
            'interpreter','latex',...
            'LineWidth', .5, ...
            'VerticalAlignment','top', 'HorizontalAlignment','left',...
            'Color','black','FontSize', fontsize_annotation,...
            'BackgroundColor','white');
        
        
        set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
        set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
        
        set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
        
        % set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
        set(gcf,'units','centimeters','OuterPosition',[0 0 2/3*columnwidth_cm .8*columnwidth_cm]);
        
        figname=sprintf('%s/%s_%s_Peaks_N_%d_q_%.3g',basedir,curmodel,varname,sqrtN,q_cur);
        fprintf('Creating figure %s\n',figname)
        if(saveswitch == 1)
           exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end

        set(gca, 'Yscale', 'log')
        ylim([0 3*max(ft_max_vals(i_N,:,ind_q))]);
        figname=sprintf('%s/%s_%s_PeaksLog_N_%d_q_%.3g',basedir,curmodel,varname,sqrtN,q_cur);
        fprintf('Creating figure %s\n',figname)
        if(saveswitch == 1)
           exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end

        set(gca, 'Yscale', 'lin')
        ylim([0 2*max(ft_sidepeak_max(i_N,:,ind_q))]);
        figname=sprintf('%s/%s_%s_PeaksZoomed_N_%d_q_%.3g',basedir,curmodel,varname,sqrtN,q_cur);
        fprintf('Creating figure %s\n',figname)
        if(saveswitch == 1)
           exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
           exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
    end
end
