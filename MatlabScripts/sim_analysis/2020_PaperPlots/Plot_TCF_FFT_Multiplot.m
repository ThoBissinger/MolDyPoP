%% 0 Initialization
clear
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/TCF/Multiplot',fig_base);

%% 0 Preparing data
sqrtN_vals = [16 32 64 128];
L_vals = [9.25 18.5 37 74];

sqrtN=128;
i_q = 6;

% T_vals = [.14 .17 .185 .20]; T_dirs = {"T_.14" "T_.17" "T_.185" "T_.20"};
T_vals = [.14 .17 .185 .20];
T_dirs = {"T_.14" "T_.17" "T_.185" "T_.20"};
N_T = numel(T_vals);
dt_sim = .01;

T_dirs_Peaks = {"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52" };
T_vals_Peaks = [.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52 ];

% T_dirs_Peaks = ["T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.15"
% "T_.17" "T_.19" "T_.21" ]; T_vals_Peaks = [.03 .05 .07 .09 .11 .13 .15
% .17 .19 .21];
basedir_mxy_Peaks=sprintf('/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_%d',sqrtN);
basedir_fmxy_Peaks=sprintf('/data/scc/thobi/210715_LinearTimeSampling/fmxy/sqrtN_%d',sqrtN);
sampfilename_Peaks='samp_Dynamics';

y_border = .008;
x_border = .008;
labelbox_size = [.025 .03];

mxy_dir=sprintf("/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_%d",sqrtN);
mxy_shorttime_dir=sprintf("/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00/sqrtN_%d",sqrtN);

fmxy_dir=sprintf("/data/scc/thobi/211201_LongerTime/fmxy/sqrtN_%d",sqrtN);
fmxy_shorttime_dir=sprintf("/data/scc/thobi/220201_ReducedSmapleStepDeltat/fmxy/sqrtN_%d",sqrtN);
% fmxy_dir_shorttime="/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00";

mxy_runmax=500;
fmxy_runmax=125;

savefile=sprintf('data_mxy_TCF_Multiplot_%d',sqrtN);
mxy_res_factor=3/4*log(mxy_runmax);
fmxy_res_factor=3/4*log(fmxy_runmax);
res_function = resolution_Gauss;
% res_function=@(t,tau) resolution_Laplace_pleateau(t,tau,tau);

marker_types=["v" "^" "o" "d"];
marker_types=["s" "s" "s" "s"];

%% The fit
mxy_mperp_cf=cell(size(T_vals));
mxy_qbin=cell(size(T_vals));
mxy_t=cell(size(T_vals));
mxy_w_cf=cell(size(T_vals));
mxy_m_cf=cell(size(T_vals));

mxy_qbin_shorttime=cell(size(T_vals));
mxy_t_shorttime=cell(size(T_vals));
mxy_mpar_cf=cell(size(T_vals));



fmxy_mperp_cf=cell(size(T_vals));
fmxy_qbin=cell(size(T_vals));
fmxy_t=cell(size(T_vals));
fmxy_w_cf=cell(size(T_vals));
fmxy_m_cf=cell(size(T_vals));
fmxy_mpar_cf=cell(size(T_vals));

fmxy_qbin_shorttime=cell(size(T_vals));
fmxy_t_shorttime=cell(size(T_vals));
fmxy_mpar_cf=cell(size(T_vals));


for i_T = 1:N_T
    T = T_vals(i_T);
    fprintf('T = %.3f\n',T);
    mxy_curfile=sprintf('%s/%s/samp_Dynamics_collect.mat',mxy_dir,T_dirs{i_T});
    mxy_curfile_st=sprintf('%s/%s/samp_Dynamics.mat',mxy_shorttime_dir,T_dirs{i_T});
    fmxy_curfile=sprintf('%s/%s/samp_Dynamics_collect.mat',fmxy_dir,T_dirs{i_T});
    fmxy_curfile_st=sprintf('%s/%s/samp_Dynamics.mat',fmxy_shorttime_dir,T_dirs{i_T});

    load(mxy_curfile,'averaging_times','qbin','gmperpmperp_collect',...
        'gww','gxx','gyy');
    mxy_mperp_cf{i_T} = real(gmperpmperp_collect)/sqrtN^2;
    mxy_w_cf{i_T} = real(gww)/sqrtN^2;
    mxy_m_cf{i_T} = (real(gxx) + real(gyy))/sqrtN^2;
    mxy_qbin{i_T} = qbin;
    mxy_t{i_T} = averaging_times;
    if (isfile(mxy_curfile_st))
        load(mxy_curfile_st,'averaging_times','qbin','gmparmpar');
    else
        load(mxy_curfile,'averaging_times','qbin','gmparmpar');
    end
    mxy_qbin_shorttime{i_T} = qbin;
    mxy_t_shorttime{i_T} = averaging_times;
    mxy_mpar_cf{i_T} = real(gmparmpar)/sqrtN^2;
    

    load(fmxy_curfile,'averaging_times','qbin','gmperpmperp_collect',...
        'gww','gxx','gyy','gmparmpar');
    fmxy_mperp_cf{i_T} = real(gmperpmperp_collect)/sqrtN^2;
    fmxy_w_cf{i_T} = real(gww)/sqrtN^2;
    fmxy_m_cf{i_T} = (real(gxx) + real(gyy))/sqrtN^2;
    fmxy_qbin{i_T} = qbin;
    fmxy_t{i_T} = averaging_times;
    if (isfile(fmxy_curfile_st))
        load(fmxy_curfile_st,'averaging_times','qbin','gmparmpar');
    else
        load(fmxy_curfile,'averaging_times','qbin','gmparmpar');
    end
    fmxy_qbin_shorttime{i_T} = qbin;
    fmxy_t_shorttime{i_T} = averaging_times;
    fmxy_mpar_cf{i_T} = real(gmparmpar)/sqrtN^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Peak Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_period = 4;   % number of periods included in DO fit
Peaks_mperp_cf=cell(2,numel(T_vals_Peaks));
Peaks_qbin=cell(2,numel(T_vals_Peaks));
Peaks_t=cell(2,numel(T_vals_Peaks));
gamma_vals_Peaks=zeros(2,numel(T_vals_Peaks));
omega_1_vals_Peak=zeros(2,numel(T_vals_Peaks));
chi_vals_Peaks=zeros(2,numel(T_vals_Peaks));
om_vals_Peaks=cell(2,numel(T_vals_Peaks));
ft_vals_Peaks=cell(2,numel(T_vals_Peaks));
om_max_vals_Peaks=zeros(2,numel(T_vals_Peaks));
ft_max_vals_Peaks=zeros(2,numel(T_vals_Peaks));
for i_T=1:numel(T_vals_Peaks)
    T=T_vals_Peaks(i_T);
    fprintf('T = %.3f\n',T);
    for i_model=1:2
        if i_model == 1
            curfile=sprintf('%s/%s/samp_Dynamics.mat',basedir_mxy_Peaks,T_dirs_Peaks{i_T});
            res_factor=mxy_res_factor;
        else
            curfile=sprintf('%s/%s/samp_Dynamics.mat',basedir_fmxy_Peaks,T_dirs_Peaks{i_T});
            res_factor=fmxy_res_factor;
        end
    
        load(curfile,'averaging_times','qbin','gmperpmperp');
        Peaks_mperp_cf{i_model,i_T} = real(gmperpmperp);
        Peaks_qbin{i_model,i_T} = qbin;
        Peaks_t{i_model,i_T} = averaging_times;
        
        t = averaging_times;
        cf = real(gmperpmperp);
        n_t=numel(t);
        n_q = numel(qbin);
        q_indices = (i_q):n_q:n_q*n_t;
        cf = cf(q_indices);
        c = fit_DampedOscillator_RealSpace(t,cf/cf(1),n_period,weightexp,'omega_1');

        gamma_vals_Peaks(i_model,i_T)=c(1);
        omega_1_vals_Peak(i_model,i_T)=c(2);
        chi_vals_Peaks(i_model,i_T)=cf(1);

        tau=res_factor/c(1);
        while tau > max(t)
            dt=(t(end)-t(1))/(numel(t)-1);
%             tau max(t)
            [t,cf] = continue_cf(t,cf/cf(1),2,max(3*dt,omega_1_vals_Peak(i_model,i_T)/8/pi));
            cf = cf * chi_vals_Peaks(i_model,i_T);
        end
        res_vals=res_function(t,tau);
    
        [ft_vals,om_vals]=FT_correlation(t, (cf .* res_vals) /cf(1), 0);
        ft_vals_Peaks{i_model,i_T}=cf(1)*real(ft_vals);
        om_vals_Peaks{i_model,i_T}=om_vals;
        [ft_max,i_max]=max(ft_vals_Peaks{i_model,i_T});
        om_max_vals_Peaks(i_model,i_T)=abs(om_vals(i_max));
        ft_max_vals_Peaks(i_model,i_T)=ft_max;

        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_period = 4;   % number of periods included in DO fit
om_skip = .001; % minimal omega separation

cf_mean_mxy=cell(1,N_T);
mxy_t_adjusted=cell(1,N_T);
gamma_vals_mxy=zeros(size(T_vals));
omega_1_vals_mxy=zeros(size(T_vals));
chi_vals_mxy=zeros(size(T_vals));
om_max_vals_mxy=zeros(size(T_vals));
ft_max_vals_mxy=zeros(size(T_vals));
ft_zero_vals_mxy=zeros(size(T_vals));
% ft_sidepeak_max_mxy=zeros(size(T_vals));
om_vals_mxy=cell(1,N_T);
ft_mean_mxy=cell(1,N_T);
ft_err_mxy=cell(1,N_T);
ft_std_mxy=cell(1,N_T);

cf_mean_fmxy=cell(1,N_T);
fmxy_t_adjusted=cell(1,N_T);
gamma_vals_fmxy=zeros(size(T_vals));
omega_1_vals_fmxy=zeros(size(T_vals));
chi_vals_fmxy=zeros(size(T_vals));
om_max_vals_fmxy=zeros(size(T_vals));
ft_max_vals_fmxy=zeros(size(T_vals));
ft_zero_vals_fmxy=zeros(size(T_vals));
% ft_sidepeak_max_fmxy=zeros(size(T_vals));
om_vals_fmxy=cell(1,N_T);
ft_mean_fmxy=cell(1,N_T);
ft_err_fmxy=cell(1,N_T);
ft_std_fmxy=cell(1,N_T);
for i_T = 1 : N_T
    for i_model = 1:2
        T=T_vals(i_T);
        if i_model == 1
            t=mxy_t{i_T};
            qbin_vals=mxy_qbin{i_T};
            cf = mxy_mperp_cf{i_T};
            res_factor=mxy_res_factor;
            runmax=mxy_runmax;
        else
            t=fmxy_t{i_T};
            qbin_vals=fmxy_qbin{i_T};
            cf=fmxy_mperp_cf{i_T};
            res_factor=fmxy_res_factor;
            runmax=fmxy_runmax;
        end
        t = t+dt_sim;
        n_t=numel(t);
        n_q = numel(qbin_vals);
        q_indices = (i_q):n_q:n_q*n_t;
        
        q = qbin_vals(i_q);

        cf = cf(:,q_indices);
        
        while 2*pi/max(t) > om_skip
            t=[t,max(t)+t];
            cf=[cf,zeros(numel(cf(:,1)),n_t)];
            n_t=numel(t);
        end
        cf_mean = mean(cf);
        c = fit_DampedOscillator_RealSpace(t,cf_mean/cf_mean(1),n_period,weightexp,'omega_1');
        
        if i_model == 1
            cf_mean_mxy{i_T}=cf_mean;
            mxy_t_adjusted{i_T}=t;
            chi_vals_mxy(i_T) = cf_mean(1);
            gamma_vals_mxy(i_T) = c(1);
            omega_1_vals_mxy(i_T) = c(2);
        else
            cf_mean_fmxy{i_T}=cf_mean;
            fmxy_t_adjusted{i_T}=t;
            chi_vals_fmxy(i_T) = cf_mean(1);
            gamma_vals_fmxy(i_T) = c(1);
            omega_1_vals_fmxy(i_T) = c(2);
        end
        
    
        tau=res_factor/c(1);
        res_vals=res_function(t,tau);
    
        [ft_vals,om_vals]=FT_correlation(t, (cf .* res_vals)', 0);
        ft_vals=real(ft_vals)';
        ft_mean=mean(real(ft_vals));
        [ft_max,i_max]=max(ft_mean);
        i_zero = find(om_vals >= 0, 1);
        runmax=numel(ft_vals(:,1));
        if i_model == 1
            ft_mean_mxy{i_T}=ft_mean;
            ft_std_mxy{i_T}=std(real(ft_vals));
            ft_err_mxy{i_T}=std(real(ft_vals))/sqrt(runmax);
            om_vals_mxy{i_T}=om_vals;
            om_max_vals_mxy(i_T) = abs(om_vals(i_max));
            ft_max_vals_mxy(i_T) = ft_max;
            ft_zero_vals_mxy(i_T) = ft_mean(i_zero);
        else
            ft_mean_fmxy{i_T}=ft_mean;
            ft_std_fmxy{i_T}=std(real(ft_vals));
            ft_err_fmxy{i_T}=std(real(ft_vals))/sqrt(runmax);
            om_vals_fmxy{i_T}=om_vals;
            om_max_vals_fmxy(i_T) = abs(om_vals(i_max));
            ft_max_vals_fmxy(i_T) = ft_max;
            ft_zero_vals_fmxy(i_T) = ft_mean(i_zero);
        end
    end
end
save(savefile);

%% 
load(savefile)
%% 0.5 Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall settings, annotatiosn etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
figure

set(gcf,'units','centimeters','OuterPosition',[0 0 pagewidth_cm 2*columnwidth_cm]);
% c_map_T=turbo(N_T+1);c_map_T=c_map(2:end,:);
% c_map_T=parula(5);
% c_map_T=autumn(5);c_map_T=c_map_T([4:-1:1],:);
% c_map_T=autumn(100);c_map_T=c_map_T([95,70,40,1],:);
c_map_T=turbo(100);c_map_T=c_map_T([65,75,85,100],:);
% c_map_T=spring(4);


%% 1 Subplot 1: mperp cf vs t for different T, MXY MODEL
subplot(2,2,1,'replace');

t_max = 5e2;
dt_skip = 8;

% c_map=turbo(N_T+1);c_map=c_map(2:end,:);
c_map=c_map_T;
for i_T = 1 : N_T
    T=T_vals(i_T);
    t=mxy_t_adjusted{i_T};
    t = t+dt_sim;
    n_t=numel(t);
    qbin_vals=mxy_qbin{i_T};
    n_q = numel(qbin_vals);
    q = qbin_vals(i_q);
    
    dispname=sprintf('$T = %.3f$', T);
%     q_indices = (i_q):n_q:n_q*n_t; cf = mxy_mperp_cf{i_T}(:,q_indices);
%     cf_mean = mean(cf); chi_vals_mxy(i_T) = cf_mean(1);
    cf_mean = cf_mean_mxy{i_T} / chi_vals_mxy(i_T);
    cf_std = std(cf / chi_vals_mxy(i_T));
    cf_err = cf_std/sqrt(mxy_runmax);
    
%     c =
%     fit_DampedOscillator_RealSpace(t,cf_mean/cf_mean(1),n_period,weightexp,'omega_1');
%     gamma_vals_mxy(i_T) = c(1); omega_1_vals_mxy(i_T) = c(2);
    c=[gamma_vals_mxy(i_T),omega_1_vals_mxy(i_T)];
    
%     cf= real(cf)/real(cf(1));
    offset= -2.5*(i_T - 1);
    
    yline(offset,'--',...
        'HandleVisibility','off',...
        'LineWidth', .4);
    hold on;
    
    
    n_skip = find((t-t(1)) > dt_skip);
    h_sim_data(i_T) = plot(t,cf_mean + offset,...
        'LineStyle', 'none', ...
        'Marker', marker_types(i_T), 'MarkerSize', 3, ...
        'LineWidth', .8, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:),...
        'MarkerFaceColor','none','MarkerEdgeColor',c_map(i_T,:), ...
        'MarkerIndices',1:n_skip:numel(t));
%         'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k', ...

    h_fit_data(i_T) = plot(t,fitfunc_DO(t,1,c) + offset,...
        'LineStyle', '-', ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color','black');

    h_leg_data(i_T) = plot(NaN*t,NaN*t,...
        'LineStyle', '-', ...
        'LineWidth', 1, ...
        'Marker', marker_types(i_T), 'MarkerSize', 4, ...
        'DisplayName', dispname, ...
        'MarkerFaceColor','none','MarkerEdgeColor',c_map(i_T,:), ...
        'Color','black');

    text(490,offset+.4,dispname,"Color",c_map(i_T,:),...
        "VerticalAlignment","bottom","HorizontalAlignment","right",...
        'FontName', 'cmr12','FontSize', fontsize_annotation,...
        "interpreter","latex");
    
end
ylim([offset-1.5 5])
xlim([0 t_max]);
% hLegend = legend('Location', 'NorthEast','interpreter','latex',...
%     'NumColumns',1);

hXLabel = xlabel('$t$','interpreter','latex');
hYLabel = ylabel('$C_{m\perp}(q,t) / \chi_{m\perp}(q)$','interpreter','latex');
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
% set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', [], ...
    'LineWidth', .5)


annotation_str = {"\textbf{MXY}",sprintf('$N = (%d)^2$',sqrtN),sprintf('$q = %.3f$',q)};
pos = get(gca,'Position');
dim=[pos(1)+x_border, ...
    pos(2)+pos(4)-y_border-2.0*labelbox_size(2), ...
    .125,.09];
h_annot=annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');
pos_box=get(h_annot,'Position');
dim=[pos_box(1),...
    pos(2)+pos(4)-pos_box(4), ... works, but weird
    pos_box(3),...
    pos_box(4)];
set(h_annot,'Position',dim);

h_text=add_subfig_label(gca,"$(a)$","ne","lin","lin",fontsize_subfiglabels);
% SE=[max(xlim), min(ylim)];
% h_text=text(SE(1),SE(2),'$(a)$',...
%     'VerticalAlignment','bottom','HorizontalAlignment','right',...
%     'interpreter','latex',...
%     'Color','black','FontSize', fontsize_ax_labels,...
%     'BackgroundColor','none');
% pos = get(gca,'Position');
% dim=[pos(1)+pos(3), ...
%     pos(2), ...
%     .001,...
%     .001];
% h_annot=annotation('textbox',dim, 'String', '$(a)$', 'FitBoxToText','on',...
%     'interpreter','latex',...
%     'EdgeColor','none',...
%     'Units','normalized',...
%     'LineWidth', .5, ...
%     'VerticalAlignment','bottom', 'HorizontalAlignment','right',...
%     'Color','black','FontSize', fontsize_ax_labels,...
%     'BackgroundColor','none');





















%% 2 Subplot 2: mperp cf FFT vs omega for different T, MXY MODEL
subplot(2,2,2,'replace');
om_max = .08;
om_skip = .001;

c_map=c_map_T;
for i_T = 1 : N_T
    T=T_vals(i_T);
    
    hold on;
    om_vals=om_vals_mxy{i_T};
    ft_mean=ft_mean_mxy{i_T};
    ft_err=ft_err_mxy{i_T};
    n_skip = find((om_vals-om_vals(1)) > om_skip,1)-1;
    
    h_sim_data(i_T) = errorbar(om_vals(1:n_skip:end),ft_mean(1:n_skip:end),ft_err(1:n_skip:end),...
        'LineStyle', 'none', ...
        'LineWidth', 1.2, ...
        'Marker','o','MarkerSize',3,...
        'CapSize',3,...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:));
end
for i_T = 1:N_T
    T = T_vals(i_T);
    om_vals=om_vals_mxy{i_T};
    ft_mean=ft_mean_mxy{i_T};
    chi_cur=chi_vals_mxy(i_T);
%     gamma_cur=gamma_vals_mxy(i_T);
%     omega_1_cur=omega_1_vals_mxy(i_T);
    fitob=fit_DampedOscillator_FourierSpace(om_vals(:),ft_mean(:),chi_cur);
    gamma_cur=fitob.gamma;
    omega_cur=fitob.omega_0;

    

    y=fitfunc_DO_reciprocal(om_vals,chi_cur,[gamma_cur,omega_cur]);
    h_sim_data_line(i_T) = plot(om_vals,y,...
        'LineStyle', '-', ...
        'LineWidth', .5, ...
        'HandleVisibility', 'off', ...
        'Color','black');
    
    dispname=sprintf('$T = %.3f$', T);
    text(.97*om_max,6.7e3-i_T*3.5e2,dispname,"Color",c_map(i_T,:),...
        "VerticalAlignment","top","HorizontalAlignment","right",...
        'FontName', 'cmr12','FontSize', fontsize_annotation,...
        "interpreter","latex");
end
% h_children=get(gca,'Children'); set(gca,'Children',[h_children(8)
% h_children(7) h_children(6) h_children(5),...
%     h_children(4) h_children(3) h_children(2) h_children(1)]);
ylim([0 7e3])
xlim([0 om_max]);

hXLabel = xlabel('$\omega$','interpreter','latex');
% hYLabel = ylabel('$S_{m\perp}(q,\omega) /
% \chi_{m\perp}(q)$','interpreter','latex');
hYLabel = ylabel('$S_{m\perp}(q,\omega)$','interpreter','latex');
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin',...
    'LineWidth', .5)


h_text=add_subfig_label(gca,"$(b)$","ne","lin","lin",fontsize_subfiglabels);
% NE=[max(xlim), max(ylim)];
% h_text=text(NE(1),NE(2),'$(b)$',...
%     'VerticalAlignment','top','HorizontalAlignment','right',...
%     'interpreter','latex',...
%     'Color','black','FontSize', fontsize_ax_labels);
% pos = get(gca,'Position');
% dim=[pos(1)+pos(3), ...
%     pos(2)+pos(4), ...
%     .001,...
%     .001];
% h_annot=annotation('textbox',dim, 'String', '$(b)$', 'FitBoxToText','on',...
%     'interpreter','latex',...
%     'EdgeColor','none',...
%     'Units','normalized',...
%     'LineWidth', .5, ...
%     'VerticalAlignment','top', 'HorizontalAlignment','right',...
%     'Color','black','FontSize', fontsize_ax_labels,...
%     'BackgroundColor','none');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax_full=gca;
pos=get(gca,'Position');
inset_pos=[pos(1)+.005,...
    pos(2)+.4*pos(4),...
    .6*pos(3),...
    .56*pos(4)];
ax_inset = axes('OuterPosition',inset_pos);
set(ax_inset, 'FontName', 'cmr12','FontSize', fontsize_axis);
hold on;
markertypes=['s' 'o'];
x=T_vals_Peaks(2:end);
y=cell(1,2);
y{1}=om_max_vals_Peaks(1,2:end);
% y{2}=ft_max_vals_Peaks(1,2:end);
y{2}=sqrt(2)*2./ft_max_vals_Peaks(1,2:end).*chi_vals_Peaks(1,2:end);
colors=["r","b"];
for i =1:2
    plot(x,y{i},...
        'LineStyle', '-', ...
        'LineWidth', .5, ...
        'Marker',markertypes(i),...
        'MarkerSize',3,...
        'HandleVisibility', 'off', ...
        'Color',colors(i));
end
xlim([0 .35]);
ylim([3e-3 2e-1]);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','log',...
    'YTick',10.^(-10:10),...
    'LineWidth', .5)
hXLabel = xlabel('$T$','interpreter','latex');
% hYLabel =
% ylabel('$S_{m\perp}(q,\omega_{\textrm{max}})$','interpreter','latex');
% hYLabel = ylabel('$S_{m\perp}^{\textrm{max}}(q)$','interpreter','latex');
% hYLabel = ylabel('$\gamma^{-1}$','interpreter','latex'); Font
set(hXLabel, 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
% set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
i_T=find(T_vals_Peaks == .23);
% text(.23,y{2}(i_T),'\ \ $S_{m\perp}^{\textrm{max}}$',...
text(.23,y{2}(i_T),'\ \ $\sqrt{2}\gamma$',...
    'VerticalAlignment','top', ...
    'HorizontalAlignment','left',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_annotation);

i_T=find(T_vals_Peaks == .14);
text(.14,y{1}(i_T),'$\omega_{\max}$',...
    'VerticalAlignment','top', ...
    'HorizontalAlignment','right',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_annotation);



















    %% 3 Subplot 3: mperp, mpar, m, w cf FFT vs omega for different T, MXY MODEL
subplot(2,2,3,'replace');
hold on;

om_max = .08;
om_skip = .001;

c_map=linspecer(4);
% c_map = c_map(5:end,:);
i_T=2;
T=T_vals(i_T);

qbin=mxy_qbin{i_T};
t=mxy_t{i_T};
n_q=numel(qbin);
n_t=numel(t);
q_indices=i_q:n_q:n_t*n_q;

qbin=mxy_qbin_shorttime{i_T};
t=mxy_t_shorttime{i_T};
n_q=numel(qbin);
n_t=numel(t);
q_indices_shorttime=i_q:n_q:n_t*n_q;

lw=1;

ts={mxy_t{i_T},mxy_t{i_T},mxy_t_shorttime{i_T},mxy_t{i_T}};
% qbins={mxy_qbin{i_T},mxy_qbin{i_T},mxy_qbin_shorttime{i_T},mxy_qbin{i_T}};
cfs={real(mxy_m_cf{i_T}(q_indices)),...
    real(mean(mxy_mperp_cf{i_T}(:,q_indices))), ...
    real(mxy_mpar_cf{i_T}(q_indices_shorttime)),...
    real(mxy_w_cf{i_T}(q_indices))};
names=["$S_{m}(q,\omega)$",...
    "$S_{m\perp}(q,\omega)$",...
    "$S_{m\parallel}(q,\omega)$",...
    "$S_{w}(q,\omega)$"];
for i = 1:4
    t=ts{i};
    cf=cfs{i};
    c = fit_DampedOscillator_RealSpace(t,cf/cf(1),n_period,weightexp,'omega_1');
    gamma=c(1);
    tau=mxy_res_factor/gamma;
    res_vals=res_function(t,tau);
    [ft_vals,om_vals]=FT_correlation(t, (cf .* res_vals)', 0);
    ft_vals=real(ft_vals);
    h_plot{i} = plot(om_vals,ft_vals,...
        'DisplayName',names(i),...
        'Marker','none',...
        'LineStyle', '-', ...
        'LineWidth', lw, ...
        'Color',c_map(i,:));
end


ylim([1e-1 1e4])
xlim([0 om_max]);

hLegend = legend('Location', 'West','interpreter','latex',...
    'NumColumns',1);

hXLabel = xlabel('$\omega$','interpreter','latex');
hYLabel = ylabel('$S(q,\omega)$','interpreter','latex');
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)


% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','log',...
    'LineWidth', .5)



% pos = get(gca,'Position');
% dim=[pos(1)+pos(3)-x_border-labelbox_size(1), ...
%     pos(2)+y_border, ...
%     labelbox_size(1),...
%     labelbox_size(2)];
% annotation('textbox',dim, 'String', 'c', 'FitBoxToText','off',...
%     'interpreter','latex',...
%     'Units','normalized',...
%     'LineWidth', .5, ...
%     'VerticalAlignment','middle', 'HorizontalAlignment','left',...
%     'Color','black','FontSize', fontsize_annotation,...
%     'BackgroundColor','white');

annotation_str = {"\textbf{MXY}",sprintf('$T = %.3f$',T)};
pos = get(gca,'Position');
dim=[pos(1)+x_border, ...
    pos(2)+pos(4)-y_border-2.5*labelbox_size(2), ...
    .125,.09];
h_annot=annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');


h_text=add_subfig_label(gca,"$(c)$","ne","lin","log",fontsize_subfiglabels);
% SE=[max(xlim), min(ylim)];
% h_text=text(SE(1),SE(2),'$(c)$',...
%     'VerticalAlignment','bottom','HorizontalAlignment','right',...
%     'interpreter','latex',...
%     'Color','black','FontSize', fontsize_ax_labels);
% pos = get(gca,'Position');
% dim=[pos(1)+pos(3), ...
%     pos(2), ...
%     .001,...
%     .001];
% h_annot=annotation('textbox',dim, 'String', '$(c)$', 'FitBoxToText','on',...
%     'interpreter','latex',...
%     'EdgeColor','none',...
%     'Units','normalized',...
%     'LineWidth', .5, ...
%     'VerticalAlignment','bottom', 'HorizontalAlignment','right',...
%     'Color','black','FontSize', fontsize_ax_labels,...
%     'BackgroundColor','none');










%% 4 Subplot 4: mperp cf FFT vs omega for different T, BOTH MODELS
subplot(2,2,4,'replace');
hold on;


om_vals=cell(1,N_T);
ft_mean=cell(1,N_T);

n_period = 4;
om_max = .06;

c_map=c_map_T;
for i_T = [3,4]
    T=T_vals(i_T);
    for i_model = 1:2
        if i_model == 1
            om_vals=om_vals_mxy{i_T};
            ft_mean=ft_mean_mxy{i_T};
            linestyle='-';
            linewidth=2;
            colorcur=c_map(i_T,:);
        else
            om_vals=om_vals_fmxy{i_T};
            ft_mean=ft_mean_fmxy{i_T};
            linestyle='-';
            linewidth=1;
            colorcur=c_map(i_T,:);
            colorcur='black';
        end
%         ft_mean=smoothen(abs(ft_mean),2,2)
    
        h_sim_data_line(i_T) = plot(om_vals,ft_mean,...
            'LineStyle', linestyle, ...
            'LineWidth', linewidth, ...
            'HandleVisibility', 'off', ...
            'Color',colorcur);

        dispname=sprintf('$T = %.3f$', T);
        text(.97*om_max,1.3e3-i_T*2.5e2,dispname,"Color",c_map(i_T,:),...
            "VerticalAlignment","top","HorizontalAlignment","right",...
            'FontName', 'cmr12','FontSize', fontsize_annotation,...
            "interpreter","latex");
    end
end
% h_children=get(gca,'Children'); set(gca,'Children',[h_children(8)
% h_children(7) h_children(6) h_children(5),...
%     h_children(4) h_children(3) h_children(2) h_children(1)]);

xlim([0 om_max]);
ylim([0 5e3]);

hXLabel = xlabel('$\omega$','interpreter','latex');
% hYLabel = ylabel('$S_{m\perp}(q,\omega) /
% \chi_{m\perp}(q)$','interpreter','latex');
hYLabel = ylabel('$S_{m\perp}(q,\omega)$','interpreter','latex');
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin',...
    'LineWidth', .5)
% ylim([1e1 Inf]);


% pos = get(gca,'Position');
% dim=[pos(1)+pos(3)-x_border-labelbox_size(1), ...
%     pos(2)+y_border, ...
%     labelbox_size(1),...
%     labelbox_size(2)];
% annotation('textbox',dim, 'String', 'd', 'FitBoxToText','off',...
%     'interpreter','latex',...
%     'Units','normalized',...
%     'LineWidth', .5, ...
%     'VerticalAlignment','middle', 'HorizontalAlignment','left',...
%     'Color','black','FontSize', fontsize_ax_labels,...
%     'BackgroundColor','white');

pos=get(gca,'Position');
arrowdim=[pos(1)+.63*pos(3),...
    pos(2)+.45*pos(4),...
    0,-.046];
h_arrow_mxy=annotation('textarrow','Position',arrowdim,'String','MXY',...
    'Units','normalized',...
    'LineWidth',.5,...
    'VerticalAlignment','bottom', 'HorizontalAlignment','center',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_annotation,...
    'HeadStyle','plain','HeadWidth',2,'HeadLength',3);

arrowdim=[pos(1)+.71*pos(3),...
    pos(2)+.85*pos(4),...
    .017,0];
h_arrow_fmxy=annotation('textarrow','Position',arrowdim,'String','DXY',...
    'Units','normalized',...
    'LineWidth',.5,...
    'VerticalAlignment','middle', 'HorizontalAlignment','right',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_annotation,...
    'HeadStyle','plain','HeadWidth',2,'HeadLength',3);

h_text=add_subfig_label(gca,"$(d)$","ne","lin","lin",fontsize_subfiglabels);
% SE=[max(xlim), min(ylim)];
% h_text=text(SE(1),SE(2),'$(d)$',...
%     'VerticalAlignment','bottom','HorizontalAlignment','right',...
%     'interpreter','latex',...
%     'Color','black','FontSize', fontsize_ax_labels);
% pos = get(gca,'Position');
% dim=[pos(1)+pos(3), ...
%     pos(2), ...
%     .001,...
%     .001];
% h_annot=annotation('textbox',dim, 'String', '$(d)$', 'FitBoxToText','on',...
%     'interpreter','latex',...
%     'EdgeColor','none',...
%     'Units','normalized',...
%     'LineWidth', .5, ...
%     'VerticalAlignment','bottom', 'HorizontalAlignment','right',...
%     'Color','black','FontSize', fontsize_ax_labels,...
%     'BackgroundColor','none');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax_full=gca;
pos=get(gca,'Position');
inset_pos=[pos(1)+.005,...
    pos(2)+.4*pos(4),...
    .6*pos(3),...
    .56*pos(4)];
ax_inset = axes('OuterPosition',inset_pos);
set(ax_inset, 'FontName', 'cmr12','FontSize', fontsize_axis);
hold on;
markertypes=['s' 'o'];
x=cell(1,2);
y=cell(1,2);
colors=["r","b"];
for i_model =1:2
    if i_model == 1
        x{i_model}=T_vals_Peaks(2:end);
%     y{i_model}=1./gamma_vals_Peaks(i_model,2:end);
        y{i_model}=ft_max_vals_Peaks(i_model,2:end);
    else
        x{i_model}=T_vals_Peaks(4:end);
        y{i_model}=ft_max_vals_Peaks(i_model,4:end);
    end
%     plot(T_vals_Peaks(2:end),ft_max_vals_Peaks(i_model,2:end),...
%     plot(T_vals_Peaks(2:end),1./gamma_vals_Peaks(i_model,2:end),...
    plot(x{i_model},y{i_model},...
        'LineStyle', '-', ...
        'LineWidth', .5, ...
        'Marker',markertypes(i_model),...
        'MarkerSize',3,...
        'HandleVisibility', 'off', ...
        'Color',colors(i_model));
end
xlim([0 .35]);
% ylim([1e0 5e4]);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','log',...
    'YTick',10.^(-10:10),...
    'LineWidth', .5)
hXLabel = xlabel('$T$','interpreter','latex');
% hYLabel =
% ylabel('$S_{m\perp}(q,\omega_{\textrm{max}})$','interpreter','latex');
hYLabel = ylabel('$S_{m\perp}^{\textrm{max}}(q)$','interpreter','latex');
% hYLabel = ylabel('$\gamma^{-1}$','interpreter','latex'); Font
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
i_T=find(T_vals_Peaks == .14);
text(.17,y{2}(i_T),'DXY',...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','left',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_annotation);

i_T=find(T_vals_Peaks == .14);
text(.14,y{1}(i_T),'MXY',...
    'VerticalAlignment','top', ...
    'HorizontalAlignment','right',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_annotation);






% dim=[pos(1),pos(2),...
%     pos(3) pos(4)];
% h_rect=annotation('rectangle',dim,'Color','magenta',...
%     'LineStyle',':','LineWidth',1.5);






















%% Saving
figname=sprintf('%s/Compare_Corrs_sqrtN_%d_q_%.3f',basedir,sqrtN,q);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end