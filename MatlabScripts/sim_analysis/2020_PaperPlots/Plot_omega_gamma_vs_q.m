%% 0 Initialization
clear all
initialization_script
saveswitch=1;
basedir=sprintf('%s/plots/omega_gamma_plots',fig_base);
i_model = 3;
if (i_model == 1)
    curmodel="mxy";
    curtitle="MXY model";
    modelname="MXY";

    sqrtN_vals = [16 32 64 128 256];
    dir="/data/scc/thobi/211201_LongerTime/mxy_3.00";
%     dir="/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00";
    dir_more_q="/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00";

    T_vals = [.11 .14 .17 .18];
    T_dirs = {"T_.11" "T_.14" "T_.17" "T_.18"};
%     T_vals = [.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52];
%     T_dirs = {"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
    
    sampfilename="samp_Dynamics";
    
    L_vals=[9.25,18.5,37,74,148];
    q_select=[1:23];
    q_select=[1 2 3 4 5 6 7 8 10 11 13 14 15 18 21 22]; % Ignores multiple counts due to binning
elseif (i_model == 2)
    
    L_vals=[16,32,64,128,256];
    
    
elseif (i_model == 3)
    curmodel="fmxy";
    curtitle="FMXY model";
    modelname="FMXY";

    sqrtN_vals = [16 32 64 128 256];
    dir="/data/scc/thobi/211201_LongerTime/fmxy";
    dir_more_q="/data/scc/thobi/210715_LinearTimeSampling/fmxy";

    T_vals = [.11 .14 .17 .18];
    T_dirs = {"T_.11" "T_.14" "T_.17" "T_.18"};
    
    sampfilename="samp_Dynamics";
    
    L_vals=[9.25,18.5,37,74,148];
    q_select=[1:23];
    q_select=[1 2 3 4 5 6 7 8 10 11 13 14 15 18 21 22]; % Ignores multiple counts due to binning
    
elseif (i_model == 4)
    L_vals=[16, 32, 64, 128];
    
end
N_N = numel(sqrtN_vals);
N_T = numel(T_vals);
N_q = numel(q_select);
markersym='o';

%% 1 Data assembly
eta_vals = zeros(1,N_T);
for i_T = 1:N_T
    T_dir=T_dirs{i_T};
    absM_vec=zeros(1,N_N);
    for i_N = 1:N_N
        sqrtN = sqrtN_vals(i_N);
        curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_more_q,sqrtN,T_dir,sampfilename);
        load(curfile,"absM_av");
        absM_vec(i_N) = absM_av;
    end
    eta_fitob = fit_eta_Magnetization_FS(absM_vec,L_vals);
    eta_vals(i_T) = eta_fitob.eta;
end
i_N = 4;
n_period = 10;
weightexp = 1;
sqrtN = sqrtN_vals(i_N);
gamma_vals=zeros(N_T,N_q);
omega_1_vals=zeros(N_T,N_q);

max_FFT_vals=zeros(N_T,N_q);
gamma_FFT_vals=zeros(N_T,N_q);
omega_2_FFT_vals=zeros(N_T,N_q);

q_vals_collect=zeros(1,N_q);
for i_T = 1:N_T
    T = T_vals(i_T);
    T_dir = T_dirs{i_T};

    curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir,sqrtN,T_dir,sampfilename);
    if (~ isfile(curfile))
        curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_more_q,sqrtN,T_dir,sampfilename);
    end
    load(curfile,"averaging_times","gmperpmperp","qbin");
    t = averaging_times;
    n_t = numel(t);
    n_q = numel(qbin);
    for i_q = 1:N_q
        if q_select(i_q) > n_q - 1
            curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_more_q,sqrtN,T_dir,sampfilename);
            load(curfile,"averaging_times","gmperpmperp","qbin");
            t = averaging_times;
            n_t = numel(t);
            n_q = numel(qbin);
        end
        q_vals_collect(i_q) = qbin(q_select(i_q));

        cf=real(gmperpmperp(q_select(i_q):n_q:end));
        cf=cf/cf(1);

        c = fit_DampedOscillator_RealSpace(t,cf,n_period,weightexp,'omega_1');
        gamma_vals(i_T,i_q) = c(1);
        omega_1_vals(i_T,i_q) = c(2);

        [ft_vals,om_vals]=FT_correlation(t, cf, 0);
        ft_vals=real(ft_vals);
        [ft_max,i_ft_max]=max(ft_vals);
        omega_2_FFT_vals(i_T,i_q)=abs(om_vals(i_ft_max));
        max_FFT_vals(i_T,i_q)=ft_max;
    end
end
omega_diff_vals = sqrt(omega_1_vals.^2 - gamma_vals.^2/4);
omega_0_vals = sqrt(omega_1_vals.^2 + gamma_vals.^2/4);

gamma_FFT_vals_a = 2./max_FFT_vals;
gamma_FFT_vals_b = omega_0_vals.^2 .* max_FFT_vals / 2;

gamma_a_vals=zeros(1,N_T);
gamma_sigma_vals=zeros(1,N_T);
omega_0_a_vals=zeros(1,N_T);
omega_0_sigma_vals=zeros(1,N_T);
omega_1_a_vals=zeros(1,N_T);
omega_1_sigma_vals=zeros(1,N_T);
omega_2_a_vals=zeros(1,N_T);
omega_2_sigma_vals=zeros(1,N_T);
fitfunc_coeffs=@(a,sigma,x) a*x.^sigma;
ind_select=1:10;
q_vec=q_vals_collect(ind_select);
for i_T = 1:N_T
    vec_cur=gamma_vals(i_T,ind_select);
    init_sigma=2;
    init_a=vec_cur(end)/q_vec(end)^init_sigma;
    fitob = fit(q_vec(:),vec_cur(:),fittype(fitfunc_coeffs),...
        'StartPoint',[init_a,init_sigma]);
    gamma_a_vals(i_T) = fitob.a;
    gamma_sigma_vals(i_T) = fitob.sigma;

    vec_cur=omega_1_vals(i_T,ind_select);
    init_sigma=1;
    init_a=vec_cur(end)/q_vec(end)^init_sigma;
    fitob = fit(q_vec(:),vec_cur(:),fittype(fitfunc_coeffs),...
        'StartPoint',[init_a,init_sigma]);
    omega_1_a_vals(i_T) = fitob.a;
    omega_1_sigma_vals(i_T) = fitob.sigma;

    vec_cur=omega_0_vals(i_T,ind_select);
    init_sigma=1;
    init_a=vec_cur(end)/q_vec(end)^init_sigma;
    fitob = fit(q_vec(:),vec_cur(:),fittype(fitfunc_coeffs),...
        'StartPoint',[init_a,init_sigma]);
    omega_0_a_vals(i_T) = fitob.a;
    omega_0_sigma_vals(i_T) = fitob.sigma;

    vec_cur=omega_2_FFT_vals(i_T,ind_select);
    init_sigma=1;
    init_a=vec_cur(end)/q_vec(end)^init_sigma;
    fitob = fit(q_vec(:),vec_cur(:),fittype(fitfunc_coeffs),...
        'StartPoint',[init_a,init_sigma]);
    omega_2_a_vals(i_T) = fitob.a;
    omega_2_sigma_vals(i_T) = fitob.sigma;
end

























%% 2 Plot gamma vs q for a few T
figure
c_map=linspecer(N_T);
for i_T = 1:N_T
    T = T_vals(i_T);
    eta = eta_vals(i_T);
    
%     semilogy(q_vals_collect,gamma_vals(i_T,:),'s',...
%         'HandleVisibility','off',...
%         'LineWidth',2,...
%         'Color',c_map(i_T,:));
    hold on;
    y_fit=fitfunc_coeffs(gamma_a_vals(i_T),gamma_sigma_vals(i_T),q_vals_collect);
%     semilogy(q_vals_collect,y_fit,'--',...
%         'HandleVisibility','off',...
%         'LineWidth',2,...
%         'Color',c_map(i_T,:));

    dispname = sprintf('$T = %.3g$, $\\gamma \\sim q^{%.2f}$',T,gamma_sigma_vals(i_T));
%     semilogy(NaN,NaN,'s--',...
%         'DisplayName',dispname,...
%         'LineWidth',2,...
%         'Color',c_map(i_T,:));

    h_fit_data(i_T) = plot(q_vals_collect,y_fit,...
        'LineStyle', '-.', ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:));

    h_sim_data(i_T) = plot(q_vals_collect,gamma_vals(i_T,:),'s',...
        'LineStyle', 'none', ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:),...
        'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k');

    h_leg_data(i_T) = plot(NaN,NaN,...
        'LineStyle', '-.', ...
        'LineWidth', 1.3, ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'DisplayName', dispname, ...
        'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k', ...
        'Color',c_map(i_T,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust axes, legend properties etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend=legend('Location','SouthEast','Interpreter','latex',...
        'NumColumns',1);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', .5,'Xscale', 'log','Yscale', 'log')

hXLabel = xlabel('$q$','interpreter','latex');
hYLabel = ylabel('$\gamma$','interpreter','latex');
h_axis = gca;

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)


annotation_str = {modelname,sprintf('$N = (%d)^2$',sqrtN)};
dim=[.685 .17+.04*N_T .1 .2];
h_annotation = annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

figname=sprintf('%s/%s_gamma_vs_q_log',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end


%% 2.1 Setting it to linear


for i_T = 1:N_T
    h_fit_data(i_T).YData = h_fit_data(i_T).YData ./ h_fit_data(i_T).XData.^2;
    h_sim_data(i_T).YData = h_sim_data(i_T).YData ./ h_sim_data(i_T).XData.^2;
end
% xlim auto;
% ylim auto;
xlim([0 .67])
ylim([0 .35]);
set(hYLabel,'String','$\gamma / q^2$');
set(hLegend,'Location','NorthEast');
set(gca, 'Xscale', 'lin','Yscale', 'lin')
dim=[.29 .725 .1 .2];
set(h_annotation,'Position',dim,'FitBoxToText','on')

figname=sprintf('%s/%s_gamma_vs_q_lin',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end































%% 3 Plot omega_1 vs q for a few T
figure
c_map=linspecer(N_T);
markersym='o';
for i_T = 1:N_T
    T = T_vals(i_T);
    eta = eta_vals(i_T);
    
    hold on;
    y_fit=fitfunc_coeffs(omega_1_a_vals(i_T),omega_1_sigma_vals(i_T),q_vals_collect);
    dispname = sprintf('$T = %.3g$, $\\omega_1 \\sim q^{%.2f}$',T,omega_1_sigma_vals(i_T));

    h_fit_data(i_T) = plot(q_vals_collect,y_fit./q_vals_collect,...
        'LineStyle', '-.', ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:));

    h_sim_data(i_T) = plot(q_vals_collect,omega_1_vals(i_T,:)./q_vals_collect,'s',...
        'LineStyle', 'none', ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:),...
        'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k');

    h_leg_data(i_T) = plot(NaN,NaN,...
        'LineStyle', '-.', ...
        'LineWidth', 1.3, ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'DisplayName', dispname, ...
        'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k', ...
        'Color',c_map(i_T,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust axes, legend properties etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend=legend('Location','NorthEast','Interpreter','latex',...
        'NumColumns',1);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')

hXLabel = xlabel('$q$','interpreter','latex');
hYLabel = ylabel('$\omega_1/q$','interpreter','latex');
h_axis = gca;
ylim([.177 .273]);

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)


annotation_str = {modelname,sprintf('$N = (%d)^2$',sqrtN)};
dim=[.16 .73 .1 .2];
h_annotation = annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

figname=sprintf('%s/%s_omega_1_vs_q_lin',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

%% 3.1 Setting it to logarithmic


for i_T = 1:N_T
    h_fit_data(i_T).YData = h_fit_data(i_T).YData .* h_fit_data(i_T).XData;
    h_sim_data(i_T).YData = h_sim_data(i_T).YData .* h_sim_data(i_T).XData;
end
% xlim auto;
% ylim auto;
xlim([0 .67])
ylim auto;
set(hYLabel,'String','$\omega_1$');
set(hLegend,'Location','SouthEast');
set(gca, 'Xscale', 'log','Yscale', 'log')
dim=[.29 .1 .1 .2];
set(h_annotation,'Position',dim,'FitBoxToText','on')

figname=sprintf('%s/%s_omega_1_vs_q_log',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end















%% 4 Plot omega_2 vs q for a few T
figure
c_map=linspecer(N_T);
markersym='o';
for i_T = 1:N_T
    T = T_vals(i_T);
    eta = eta_vals(i_T);
    
    hold on;
    y_fit=fitfunc_coeffs(omega_2_a_vals(i_T),omega_2_sigma_vals(i_T),q_vals_collect);
    dispname = sprintf('$T = %.3g$, $\\omega_2 \\sim q^{%.2f}$',T,omega_2_sigma_vals(i_T));

    h_fit_data(i_T) = plot(q_vals_collect,y_fit./q_vals_collect,...
        'LineStyle', '-.', ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:));

    h_sim_data(i_T) = plot(q_vals_collect,omega_2_FFT_vals(i_T,:)./q_vals_collect,'s',...
        'LineStyle', 'none', ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:),...
        'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k');

    h_leg_data(i_T) = plot(NaN,NaN,...
        'LineStyle', '-.', ...
        'LineWidth', 1.3, ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'DisplayName', dispname, ...
        'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k', ...
        'Color',c_map(i_T,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust axes, legend properties etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend=legend('Location','NorthEast','Interpreter','latex',...
        'NumColumns',1);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')

hXLabel = xlabel('$q$','interpreter','latex');
hYLabel = ylabel('$\omega_2/q$','interpreter','latex');
h_axis = gca;
ylim([.177 .273]);

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)


annotation_str = {modelname,sprintf('$N = (%d)^2$',sqrtN)};
dim=[.16 .73 .1 .2];
h_annotation = annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

figname=sprintf('%s/%s_omega_2_vs_q_lin',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

%% 4.1 Setting it to logarithmic


for i_T = 1:N_T
    h_fit_data(i_T).YData = h_fit_data(i_T).YData .* h_fit_data(i_T).XData;
    h_sim_data(i_T).YData = h_sim_data(i_T).YData .* h_sim_data(i_T).XData;
end
% xlim auto;
% ylim auto;
xlim([0 .67])
ylim auto;
set(hYLabel,'String','$\omega_2$');
set(hLegend,'Location','SouthEast');
set(gca, 'Xscale', 'log','Yscale', 'log')
dim=[.29 .1 .1 .2];
set(h_annotation,'Position',dim,'FitBoxToText','on')

figname=sprintf('%s/%s_omega_2_vs_q_log',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end












%% 4 Plot omega_0 vs q for a few T
figure
c_map=linspecer(N_T);
markersym='o';
for i_T = 1:N_T
    T = T_vals(i_T);
    eta = eta_vals(i_T);
    
    hold on;
    y_fit=fitfunc_coeffs(omega_0_a_vals(i_T),omega_0_sigma_vals(i_T),q_vals_collect);
    dispname = sprintf('$T = %.3g$, $\\omega_0 \\sim q^{%.2f}$',T,omega_0_sigma_vals(i_T));

    h_fit_data(i_T) = plot(q_vals_collect,y_fit./q_vals_collect,...
        'LineStyle', '-.', ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:));

    h_sim_data(i_T) = plot(q_vals_collect,omega_0_vals(i_T,:)./q_vals_collect,'s',...
        'LineStyle', 'none', ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_T,:),...
        'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k');

    h_leg_data(i_T) = plot(NaN,NaN,...
        'LineStyle', '-.', ...
        'LineWidth', 1.3, ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'DisplayName', dispname, ...
        'MarkerFaceColor',c_map(i_T,:),'MarkerEdgeColor','k', ...
        'Color',c_map(i_T,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust axes, legend properties etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend=legend('Location','NorthEast','Interpreter','latex',...
        'NumColumns',1);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')

hXLabel = xlabel('$q$','interpreter','latex');
hYLabel = ylabel('$\omega_0/q$','interpreter','latex');
h_axis = gca;

xlim([0 .67])
% ylim([.177 .273]);

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)


annotation_str = {modelname,sprintf('$N = (%d)^2$',sqrtN)};
dim=[.29 .73 .1 .2];
h_annotation = annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

figname=sprintf('%s/%s_omega_0_vs_q_lin',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

%% 4.1 Setting it to logarithmic


for i_T = 1:N_T
    h_fit_data(i_T).YData = h_fit_data(i_T).YData .* h_fit_data(i_T).XData;
    h_sim_data(i_T).YData = h_sim_data(i_T).YData .* h_sim_data(i_T).XData;
end
% xlim auto;
% ylim auto;
xlim([0 .67])
ylim auto;
set(hYLabel,'String','$\omega_0$');
set(hLegend,'Location','NorthWest');
set(gca, 'Xscale', 'log','Yscale', 'log')
dim=[.56 .73 .1 .2];
set(h_annotation,'Position',dim,'FitBoxToText','on')

figname=sprintf('%s/%s_omega_0_vs_q_log',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

%% 4 Plot omega_0 vs q for a few T
figure
c_map=linspecer(N_T);
for i_T = 1:N_T
    T = T_vals(i_T);
    eta = eta_vals(i_T);
    
    semilogy(q_vals_collect,omega_0_vals(i_T,:),'s',...
        'HandleVisibility','off',...
        'LineWidth',2,...
        'Color',c_map(i_T,:));
    hold on;
    y_fit=fitfunc_coeffs(omega_0_a_vals(i_T),omega_0_sigma_vals(i_T),q_vals_collect);
    semilogy(q_vals_collect,y_fit,'--',...
        'HandleVisibility','off',...
        'LineWidth',2,...
        'Color',c_map(i_T,:));

    dispname = sprintf('$T = %.3g$, $\\omega_0 \\sim q^{%.2f}$',T,omega_0_sigma_vals(i_T));
    semilogy(NaN,NaN,'s--',...
        'DisplayName',dispname,...
        'LineWidth',2,...
        'Color',c_map(i_T,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust axes, legend properties etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend=legend('Location','SouthEast','Interpreter','latex',...
        'NumColumns',1);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', .5,'Xscale', 'log','Yscale', 'lin')

hXLabel = xlabel('$q$','interpreter','latex');
hYLabel = ylabel('$\omega_0$','interpreter','latex');
h_axis = gca;

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)


annotation_str = {modelname,sprintf('$N = (%d)^2$',sqrtN)};
dim=[.685 .17+.04*N_T .1 .2];
annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

figname=sprintf('%s/%s_omega_0_vs_q',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end













return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0 Initialization (the old way)
clear all
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;


mxydata_name='mxy/rho_3.00_dynamics_LinearTime.mat';
xydata_name='xy/xy_dynamics_LinearTime.mat';
fmxydata_name='fmxy/fmxy_dynamics_LinearTime.mat';

mxyfit_name='mxy/rho_3.00_CritExpFit.mat';
xyfit_name='xy/xy_CritExpFit.mat';
fmxyfit_name='fmxy/fmxy_CritExpFit.mat';

mxyfit_TCF_q_name='mxy/rho_3.00_TCF_q_fit_LinearTime.mat';
xyfit_TCF_q_name='xy/xy_TCF_q_fit_LinearTime.mat';
fmxyfit_TCF_q_name='fmxy/fmxy_TCF_q_fit_LinearTime.mat';

mxyfit_FS_name='mxy/rho_3.00_DataCollapse_SCF.mat';
xyfit_FS_name='xy/xy_DataCollapse_SCF.mat';
fmxyfit_FS_name='fmxy/fmxy_DataCollapse_SCF.mat';

for i_model = [1]

    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=load(mxydata_name);
        fitdata=load(mxyfit_TCF_q_name);
        critexpfit=load(mxyfit_name);
        
        L_vals=[9.25,18.5,37,74,148];
%         T_select=[9:2:13,14:21];
%         T_select=1:36;
        T_select_gamma=[3:3:18];
        T_select_omega=[1:3:12,13:15];

        T_max = .22;

    elseif (i_model == 2)
        curmodel="xy";
        curtitle="SXY model";
        
        data=load(xydata_name);
        fitdata=load(xyfit_TCF_q_name);
        critexpfit=load(xyfit_name);

%         T_select=[17:3:38];
        T_select_gamma=[3:3:18];
        T_select_omega=[1:3:12,13:15];

        L_vals=[16,32,64,128,256];
        
        
    elseif (i_model == 3)
        curmodel="fmxy";
        curtitle="FMXY model";
        
        data=load(fmxydata_name);
        fitdata=load(fmxyfit_TCF_q_name);
        critexpfit=load(fmxyfit_name);

        L_vals=[9.25,18.5,37,74,148];
%         T_select=1:36;
        T_select_gamma=[3:3:18];
        T_select_omega=[1:3:12,13:15];

    end

    T_vals=data.('T_vals');
    sqrtN_vals=data.('sqrtN_vals');
%     gmperpmperp=data.('gmperpmperp');
    TCF_times=data.('TCF_times');
    qbin = fitdata.('qbin');
    
    param_TCFSpin_omega_1_q_DO = fitdata.('param_TCFSpin_omega_1_q_DO');
    omega_1_a = fitdata.('omega_1_a');
    omega_1_sigma = fitdata.('omega_1_sigma');
    gamma_a = fitdata.('gamma_a');
    gamma_sigma = fitdata.('gamma_sigma');
    gamma_a_with_offset = fitdata.('gamma_a_with_offset');
    gamma_sigma_with_offset = fitdata.('gamma_sigma_with_offset');
    gamma_c_with_offset = fitdata.('gamma_c_with_offset');
    
    
    eta_vals = critexpfit.('eta_vals');
    fit_T = critexpfit.('T_vals');

    fitfunc_DO=@(times,corrfunc_1,c) corrfunc_1 * exp(-c(1) * times/2) .* (cos(c(2) * times) + .5*c(1)/c(2)*sin(c(2) * times));  
    
    i_N = numel(sqrtN_vals)-1;
    q_vals = qbin{i_N,1};


%     [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% gamma plot
    T_select = T_select_gamma;
    N_sqrtN = numel(sqrtN_vals);
    c_map = linspecer(numel(T_select));
    
    figure(100*i_model + 1)
    for ind_T = 1:length(T_select)
        
        i_T = T_select(ind_T);
        T = T_vals(i_T);
        
        dispname=sprintf('T = %.3f', T);    
        
        curcoeff = param_TCFSpin_omega_1_q_DO{i_N,i_T};
        curgamma = curcoeff(1:2:end);
        curomega = curcoeff(2:2:end);
        
        plot(q_vals,curgamma,...
            'Color',c_map(ind_T,:),'DisplayName',dispname,...
            'LineStyle', '--', 'Marker','o', ...
            'LineWidth',1.8);
        hold on;
    end

    h_axis = gca;
    set(h_axis,'xscale','log');

    hXLabel = xlabel('$q$','interpreter','latex');
    hYLabel = ylabel('$\gamma(q)$','interpreter','latex');

%     xlim([q_vals(2),Inf]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Legend, axes etc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hLegend = legend('Location', 'Best','interpreter','latex','FontSize', 10,...
        'NumColumns',2);

    hTitle = title(curtitle);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Font
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set(h_axis, 'FontSize', 12)
    set([hXLabel, hYLabel], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Adjust axes properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid','off',...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'log','Yscale', 'log')
    xlim([.8*q_vals(2),Inf]);

    figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/omega_gamma_plots/%s_SpinTCF_sqrtN_%d_gamma',curmodel,sqrtN_vals(i_N));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    
    set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid','off',...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'lin','Yscale', 'log')
    xlim([0,Inf]);
    legend('Location', 'Best');

    figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/omega_gamma_plots/%s_SpinTCF_sqrtN_%d_gamma_linlog',curmodel,sqrtN_vals(i_N));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid','off',...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')
    legend('Location', 'Best');
    ylim([-Inf,.005]);

    figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/omega_gamma_plots/%s_SpinTCF_sqrtN_%d_gamma_linlin',curmodel,sqrtN_vals(i_N));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    
    %% Gamma/q
    figure(100*i_model + 2)
    for ind_T = 1:length(T_select)
        
        i_T = T_select(ind_T);
        T = T_vals(i_T);
        
        dispname=sprintf('T = %.3f', T);    
        
        curcoeff = param_TCFSpin_omega_1_q_DO{i_N,i_T};
        curgamma = curcoeff(1:2:end);
        curomega = curcoeff(2:2:end);
        
        plot(q_vals,curgamma./q_vals,...
            'Color',c_map(ind_T,:),'DisplayName',dispname,...
            'LineStyle', '--', 'Marker','o', ...
            'LineWidth',1.8);
        hold on;
    end

    h_axis = gca;
    set(h_axis,'xscale','log');

    hXLabel = xlabel('$q$','interpreter','latex');
    hYLabel = ylabel('$\gamma(q)/q$','interpreter','latex');

%     xlim([q_vals(2),Inf]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Legend, axes etc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hLegend = legend('Location', 'Best','interpreter','latex','FontSize', 10,...
        'NumColumns',2);

    hTitle = title(curtitle);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Font
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set(h_axis, 'FontSize', 12)
    set([hXLabel, hYLabel], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Adjust axes properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid','off',...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'log','Yscale', 'log')


    figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/omega_gamma_plots/%s_SpinTCF_sqrtN_%d_gamma_by_q',curmodel,sqrtN_vals(i_N));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    
    
    
    
    %% gamma plot with fit
    T_select = T_select_gamma;
    figure(100*i_model + 3)
    hold on;
    for ind_T = 1:length(T_select)
        
        i_T = T_select(ind_T);
        T = T_vals(i_T);
        
        dispname=sprintf('Fit, T = %.3f', T);    
        
        curcoeff = param_TCFSpin_omega_1_q_DO{i_N,i_T};
        curgamma = curcoeff(1:2:end);
        curomega = curcoeff(2:2:end);
        
        plot(q_vals,abs(curgamma),...
            'Color',c_map(ind_T,:),'HandleVisibility','off',...
            'LineStyle', 'none', 'Marker','o', 'MarkerSize',8,...
            'LineWidth',1.8);
%         plot(q_vals,gamma_a_with_offset(i_N,i_T)*q_vals.^(gamma_sigma_with_offset(i_N,i_T)+gamma_c_with_offset(i_N,i_T)),...
        plot(q_vals,gamma_a(i_N,i_T)*q_vals.^(gamma_sigma(i_N,i_T)),...
            'Color',c_map(ind_T,:),'DisplayName',dispname,...
            'LineStyle', '--', ...
            'LineWidth',1.8);
    end

    h_axis = gca;
    set(h_axis,'xscale','log');

    hXLabel = xlabel('$q$','interpreter','latex');
    hYLabel = ylabel('$\gamma(q)$','interpreter','latex');

%     ylim([3e-3,Inf]);
%     xlim([q_vals(2),Inf]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Legend, axes etc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hLegend = legend('Location', 'Best','interpreter','latex','FontSize', 10,...
        'NumColumns',2);

    hTitle = title(curtitle);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Font
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set(h_axis, 'FontSize', 12)
    set([hXLabel, hYLabel], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Adjust axes properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid','off',...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'log','Yscale', 'log')


    figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/omega_gamma_plots/%s_SpinTCF_sqrtN_%d_gamma_fit',curmodel,sqrtN_vals(i_N));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    
    
    

    %% omega plot
    T_select = T_select_omega;
    c_map = linspecer(numel(T_select));
    figure(100*i_model + 11)
    hold on;
    for ind_T = 1:length(T_select)
        
        i_T = T_select(ind_T);
        T = T_vals(i_T);
        
        dispname=sprintf('T = %.3f', T);    
        
        curcoeff = param_TCFSpin_omega_1_q_DO{i_N,i_T};
        curgamma = curcoeff(1:2:end);
        curomega = curcoeff(2:2:end);
        
        plot(q_vals,abs(curomega),...
            'Color',c_map(ind_T,:),'DisplayName',dispname,...
            'LineStyle', '--', 'Marker','o', ...
            'LineWidth',1.8);
    end

    h_axis = gca;
    set(h_axis,'xscale','log');

    hXLabel = xlabel('$q$','interpreter','latex');
    hYLabel = ylabel('$\omega(q)$','interpreter','latex');

%     ylim([3e-3,Inf]);
%     xlim([q_vals(2),Inf]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Legend, axes etc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hLegend = legend('Location', 'Best','interpreter','latex','FontSize', 10,...
        'NumColumns',2);

    hTitle = title(curtitle);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Font
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set(h_axis, 'FontSize', 12)
    set([hXLabel, hYLabel], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Adjust axes properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid','off',...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'log','Yscale', 'log')
    xlim([.8*q_vals(2),Inf]);


    figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/omega_gamma_plots/%s_SpinTCF_sqrtN_%d_omega',curmodel,sqrtN_vals(i_N));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    
    set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid','off',...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'lin','Yscale', 'log')
    xlim([0,Inf]);
    hLegend = legend('Location', 'Best','interpreter','latex','FontSize', 10,...
        'NumColumns',2);


    figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/omega_gamma_plots/%s_SpinTCF_sqrtN_%d_omega_linlog',curmodel,sqrtN_vals(i_N));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    
    set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid','off',...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')
    hLegend = legend('Location', 'Best','interpreter','latex','FontSize', 10,...
        'NumColumns',2);

    figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/omega_gamma_plots/%s_SpinTCF_sqrtN_%d_omega_linlin',curmodel,sqrtN_vals(i_N));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    %% omega/q^(1-eta/2) plot
    
    figure(100*i_model + 12)
    T_select = T_select_omega;
    
    for ind_T = 1:length(T_select)
        
        i_T = T_select(ind_T);
        T = T_vals(i_T);
        
        dispname=sprintf('T = %.3f', T);    
        
        curcoeff = param_TCFSpin_omega_1_q_DO{i_N,i_T};
        curgamma = curcoeff(1:2:end);
        curomega = curcoeff(2:2:end);
        
        cureta = eta_vals(end,find(fit_T == T,1));
        
        cureta = eta_vals(end,find(fit_T == T,1));
        
        plot(q_vals,abs(curomega)./q_vals.^(1-cureta/2),...
            'Color',c_map(ind_T,:),'DisplayName',dispname,...
            'LineStyle', '--', 'Marker','o', ...
            'LineWidth',1.8);
        hold on;
    end

    h_axis = gca;
    set(h_axis,'xscale','log');

    hXLabel = xlabel('$q$','interpreter','latex');
    hYLabel = ylabel('$\omega(q)/q^{1-\eta/2}$','interpreter','latex');

%     ylim([3e-3,Inf]);
%     xlim([q_vals(2),Inf]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Legend, axes etc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hLegend = legend('Location', 'Best','interpreter','latex','FontSize', 10,...
        'NumColumns',2);

    hTitle = title(curtitle);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Font
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set(h_axis, 'FontSize', 12)
    set([hXLabel, hYLabel], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Adjust axes properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid','off',...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'log','Yscale', 'log')


    figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/omega_gamma_plots/%s_SpinTCF_sqrtN_%d_omega_by_qto1minuseta',curmodel,sqrtN_vals(i_N));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    

    %% omega/q plot
    
    figure(100*i_model + 13)
    T_select = T_select_omega;
    
    for ind_T = 1:length(T_select)
        
        i_T = T_select(ind_T);
        T = T_vals(i_T);
        
        dispname=sprintf('T = %.3f', T);    
        
        curcoeff = param_TCFSpin_omega_1_q_DO{i_N,i_T};
        curgamma = curcoeff(1:2:end);
        curomega = curcoeff(2:2:end);
        
        cureta = eta_vals(end,find(fit_T == T,1));
        
        cureta = eta_vals(end,find(fit_T == T,1));
        
        plot(q_vals,abs(curomega)./q_vals,...
            'Color',c_map(ind_T,:),'DisplayName',dispname,...
            'LineStyle', '--', 'Marker','o', ...
            'LineWidth',1.8);
        hold on;
    end

    h_axis = gca;
    set(h_axis,'xscale','log');

    hXLabel = xlabel('$q$','interpreter','latex');
    hYLabel = ylabel('$\omega(q)/q$','interpreter','latex');

%     ylim([3e-3,Inf]);
%     xlim([q_vals(2),Inf]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Legend, axes etc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hLegend = legend('Location', 'Best','interpreter','latex','FontSize', 10,...
        'NumColumns',2);

    hTitle = title(curtitle);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Font
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set(h_axis, 'FontSize', 12)
    set([hXLabel, hYLabel], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Adjust axes properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid','off',...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'log','Yscale', 'log')


    figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/omega_gamma_plots/%s_SpinTCF_sqrtN_%d_omega_by_q',curmodel,sqrtN_vals(i_N));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    
    
    
    
    %% omega plot with fit
    T_select = T_select_omega;
    figure(100*i_model + 14)
    hold on;
    for ind_T = 1:length(T_select)
        
        i_T = T_select(ind_T);
        T = T_vals(i_T);
        
        dispname=sprintf('Fit, T = %.3f', T);    
        
        curcoeff = param_TCFSpin_omega_1_q_DO{i_N,i_T};
        curgamma = curcoeff(1:2:end);
        curomega = curcoeff(2:2:end);
        
        plot(q_vals,abs(curomega),...
            'Color',c_map(ind_T,:),'HandleVisibility','off',...
            'LineStyle', 'none', 'Marker','o', 'MarkerSize',8,...
            'LineWidth',1.8);
        plot(q_vals,omega_1_a(i_N,i_T)*q_vals.^(omega_1_sigma(i_N,i_T)),...
            'Color',c_map(ind_T,:),'DisplayName',dispname,...
            'LineStyle', '--', ...
            'LineWidth',1.8);
    end

    h_axis = gca;
    set(h_axis,'xscale','log');

    hXLabel = xlabel('$q$','interpreter','latex');
    hYLabel = ylabel('$\omega(q)$','interpreter','latex');

%     ylim([3e-3,Inf]);
    
%     xlim([q_vals(2),Inf]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Legend, axes etc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hLegend = legend('Location', 'Best','interpreter','latex','FontSize', 10,...
        'NumColumns',2);

    hTitle = title(curtitle);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Font
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set(h_axis, 'FontSize', 12)
    set([hXLabel, hYLabel], 'FontSize', 12)
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Adjust axes properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h_axis, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid','off',...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'log','Yscale', 'log')


    figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/omega_gamma_plots/%s_SpinTCF_sqrtN_%d_omega_fit',curmodel,sqrtN_vals(i_N));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end
    