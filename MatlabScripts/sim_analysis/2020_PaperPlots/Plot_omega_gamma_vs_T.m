%% 0 Initialization
clear all
initialization_script
saveswitch=1;
basedir=sprintf('%s/plots/omega_gamma_plots',fig_base);
i_model = 1;
if (i_model == 1)
    curmodel="mxy";
    curtitle="MXY model";
    modelname="MXY";

    sqrtN_vals = [16 32 64 128];
    dir="/data/scc/thobi/211201_LongerTime/mxy_3.00";
    dir_more_q="/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00";

    T_vals = [.03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25];
    T_dirs = {"T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25"};
%     T_vals = [.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52];
%     T_dirs = {"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
    
    sampfilename="samp_Dynamics";
    
    L_vals=[9.25,18.5,37,74,148];
%     q_select=[1,3,6,9,12];
    q_select=[1,3,6,9];

elseif (i_model == 2)
    
    L_vals=[16,32,64,128,256];
    
    
elseif (i_model == 3)
    curmodel="fmxy";
    curtitle="FMXY model";
    modelname="FMXY";

    sqrtN_vals = [16 32 64 128 256];
    dir="/data/scc/thobi/211201_LongerTime/fmxy";
    dir_more_q="/data/scc/thobi/210715_LinearTimeSampling/fmxy";

    T_vals = [.03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25];
    T_dirs = {"T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25"};
    
    sampfilename="samp_Dynamics";
    
    L_vals=[9.25,18.5,37,74,148];
    q_select=[1,3,6,9];
    
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
    eta_fitob = fit_eta_Magnetization_FS(absM_vec(1:N_N),L_vals(1:N_N));
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
% omega_0_vals=zeros(N_T,N_q);
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
%         omega_0_vals(i_T,i_q) = sqrt(c(2)^2 - gamma_cur^2/4);
    end
end
omega_0_vals = sqrt(omega_1_vals.^2 + gamma_vals.^2/4);
gamma_FFT_vals_a = 2./max_FFT_vals;
gamma_FFT_vals_b = omega_0_vals.^2 .* max_FFT_vals / 2;

omega_0_combi_vals = real(sqrt(2 * omega_1_vals.^2 - omega_2_FFT_vals.^2));


gamma_a_vals=zeros(1,N_T);
gamma_sigma_vals=zeros(1,N_T);
omega_0_a_vals=zeros(1,N_T);
omega_0_sigma_vals=zeros(1,N_T);
omega_0_combi_a_vals=zeros(1,N_T);
omega_0_combi_sigma_vals=zeros(1,N_T);
omega_1_a_vals=zeros(1,N_T);
omega_1_sigma_vals=zeros(1,N_T);
c_spinwavespeed_vals=zeros(1,N_T);
fitfunc_coeffs=@(a,sigma,x) a*x.^sigma;
fitfunc_spinwave=@(c,x) c*x;
ind_select=1:numel(q_vals_collect);
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

    vec_cur=omega_1_vals(i_T,ind_select);
    init_c = .2;
    fitob=fit(q_vec(:),vec_cur(:),fittype(fitfunc_spinwave),...
        'StartPoint',[init_c]);
    c_spinwavespeed_vals(i_T) = fitob.c;

    vec_cur=omega_0_vals(i_T,ind_select);
    init_sigma=1;
    init_a=vec_cur(end)/q_vec(end)^init_sigma;
    fitob = fit(q_vec(:),vec_cur(:),fittype(fitfunc_coeffs),...
        'StartPoint',[init_a,init_sigma]);
    omega_0_a_vals(i_T) = fitob.a;
    omega_0_sigma_vals(i_T) = fitob.sigma;

    vec_cur=omega_0_combi_vals(i_T,ind_select);
    init_sigma=1;
    init_a=vec_cur(end)/q_vec(end)^init_sigma;
    fitob = fit(q_vec(:),vec_cur(:),fittype(fitfunc_coeffs),...
        'StartPoint',[init_a,init_sigma]);
    omega_0_combi_a_vals(i_T) = fitob.a;
    omega_0_combi_sigma_vals(i_T) = fitob.sigma;
end


























%% 2 Plot gamma/q^2 vs T for a few q
figure
c_map=linspecer(N_q);
for i_q = 1:N_q
    q = q_vals_collect(i_q);
    dispname = sprintf('$q = %.2f$',q);
%     semilogy(T_vals,gamma_vals(:,i_q)/q^2,'s--',...
%         'DisplayName',dispname,...
%         'LineWidth',2,...
%         'Color',c_map(i_q,:));
    hold on;
    data_cur = gamma_vals(:,i_q);
    h_fit_data(i_q) = plot(T_vals,data_cur/q^2,...
        'LineStyle', '-.', ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_q,:));

    h_sim_data(i_q) = plot(T_vals,data_cur/q^2,'s',...
        'LineStyle', 'none', ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_q,:),...
        'MarkerFaceColor',c_map(i_q,:),'MarkerEdgeColor','k');

    h_leg_data(i_q) = plot(NaN,NaN,...
        'LineStyle', '-.', ...
        'LineWidth', 1.3, ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'DisplayName', dispname, ...
        'MarkerFaceColor',c_map(i_q,:),'MarkerEdgeColor','k', ...
        'Color',c_map(i_q,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust axes, legend properties etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend=legend('Location','NorthWest','Interpreter','latex',...
        'NumColumns',1);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', .5,'Xscale', 'lin','Yscale', 'log')

hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$\gamma/q^2$','interpreter','latex');
h_axis = gca;

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)


annotation_str = {modelname,sprintf('$N = (%d)^2$',sqrtN)};
dim=[.685 .12 .1 .2];
annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

figname=sprintf('%s/%s_gamma_byq2_vs_T',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end






%% 2.1 Setting it to gamma vs T

for i_q = 1:N_q
    q = q_vals_collect(i_q);
    h_fit_data(i_q).YData = h_fit_data(i_q).YData .* q.^2;
    h_sim_data(i_q).YData = h_sim_data(i_q).YData .* q.^2;
end
set(hYLabel,'String','$\gamma$');
set(hLegend,'Location','NorthWest');
% set(gca, 'Xscale', 'lin','Yscale', 'lin')
% dim=[.29 .725 .1 .2];
% set(h_annotation,'Position',dim,'FitBoxToText','on')

figname=sprintf('%s/%s_gamma_vs_T',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end



















%% 3 Plot omega_1/q vs T for a few q
figure
c_map=linspecer(N_q);
for i_q = 1:N_q
    q = q_vals_collect(i_q);
    dispname = sprintf('$q = %.2f$',q);
%     semilogy(T_vals,gamma_vals(:,i_q)/q^2,'s--',...
%         'DisplayName',dispname,...
%         'LineWidth',2,...
%         'Color',c_map(i_q,:));
    hold on;
    data_cur = omega_1_vals(:,i_q);
    h_fit_data(i_q) = plot(T_vals,data_cur/q,...
        'LineStyle', '-.', ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_q,:));

    h_sim_data(i_q) = plot(T_vals,data_cur/q,'s',...
        'LineStyle', 'none', ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_q,:),...
        'MarkerFaceColor',c_map(i_q,:),'MarkerEdgeColor','k');

    h_leg_data(i_q) = plot(NaN,NaN,...
        'LineStyle', '-.', ...
        'LineWidth', 1.3, ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'DisplayName', dispname, ...
        'MarkerFaceColor',c_map(i_q,:),'MarkerEdgeColor','k', ...
        'Color',c_map(i_q,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust axes, legend properties etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend=legend('Location','SouthWest','Interpreter','latex',...
        'NumColumns',1);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')

hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$\omega_1/q$','interpreter','latex');
h_axis = gca;

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)


annotation_str = {modelname,sprintf('$N = (%d)^2$',sqrtN)};
dim=[.165 .17+.04*N_q .1 .2];
h_annotation=annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

figname=sprintf('%s/%s_omega_1_byq_vs_T',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end


%% 3.1 Setting it to omega_1 vs T

for i_q = 1:N_q
    q = q_vals_collect(i_q);
    h_fit_data(i_q).YData = h_fit_data(i_q).YData .* q;
    h_sim_data(i_q).YData = h_sim_data(i_q).YData .* q;
end
set(hYLabel,'String','$\omega_1$');
set(hLegend,'Location','SouthWest');
% set(gca, 'Xscale', 'lin','Yscale', 'lin')
dim=[.43 .11 .1 .2];
set(h_annotation,'Position',dim,'FitBoxToText','on')

figname=sprintf('%s/%s_omega_1_vs_T',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end


















%% 4 Plot omega_2/q vs T for a few q
figure
c_map=linspecer(N_q);
for i_q = 1:N_q
    q = q_vals_collect(i_q);
    dispname = sprintf('$q = %.2f$',q);

    hold on;
    data_cur = omega_2_FFT_vals(:,i_q);
    h_fit_data(i_q) = plot(T_vals,data_cur/q,...
        'LineStyle', '-.', ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_q,:));

    h_sim_data(i_q) = plot(T_vals,data_cur/q,'s',...
        'LineStyle', 'none', ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_q,:),...
        'MarkerFaceColor',c_map(i_q,:),'MarkerEdgeColor','k');

    h_leg_data(i_q) = plot(NaN,NaN,...
        'LineStyle', '-.', ...
        'LineWidth', 1.3, ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'DisplayName', dispname, ...
        'MarkerFaceColor',c_map(i_q,:),'MarkerEdgeColor','k', ...
        'Color',c_map(i_q,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust axes, legend properties etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend=legend('Location','SouthWest','Interpreter','latex',...
        'NumColumns',1);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')

hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$\omega_2/q$','interpreter','latex');
h_axis = gca;

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)


annotation_str = {modelname,sprintf('$N = (%d)^2$',sqrtN)};
dim=[.165 .17+.04*N_q .1 .2];
h_annotation=annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

figname=sprintf('%s/%s_omega_2_byq_vs_T',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end


%% 4.1 Setting it to omega_2 vs T

for i_q = 1:N_q
    q = q_vals_collect(i_q);
    h_fit_data(i_q).YData = h_fit_data(i_q).YData .* q;
    h_sim_data(i_q).YData = h_sim_data(i_q).YData .* q;
end
set(hYLabel,'String','$\omega_2$');
set(hLegend,'Location','SouthWest');
% set(gca, 'Xscale', 'lin','Yscale', 'lin')
dim=[.43 .11 .1 .2];
set(h_annotation,'Position',dim,'FitBoxToText','on')

figname=sprintf('%s/%s_omega_2_vs_T',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end




















%% 5 Plot omega_0/q vs T for a few q
figure
c_map=linspecer(N_q);
for i_q = 1:N_q
    q = q_vals_collect(i_q);
    dispname = sprintf('$q = %.2f$',q);
    hold on;
    data_cur = omega_0_combi_vals(:,i_q);
    h_fit_data(i_q) = plot(T_vals,data_cur/q,...
        'LineStyle', '-.', ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_q,:));

    h_sim_data(i_q) = plot(T_vals,data_cur/q,'s',...
        'LineStyle', 'none', ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_q,:),...
        'MarkerFaceColor',c_map(i_q,:),'MarkerEdgeColor','k');

    h_leg_data(i_q) = plot(NaN,NaN,...
        'LineStyle', '-.', ...
        'LineWidth', 1.3, ...
        'Marker', markersym, 'MarkerSize', 3.8, ...
        'DisplayName', dispname, ...
        'MarkerFaceColor',c_map(i_q,:),'MarkerEdgeColor','k', ...
        'Color',c_map(i_q,:));
end
% for i_q = 1:N_q
%     q = q_vals_collect(i_q);
%     dispname = sprintf('$q = %.2f$',q);
%     plot(T_vals,omega_1_vals(:,i_q)/q,'s--',...
%         'DisplayName',dispname,...
%         'LineWidth',2,...
%         'Color',c_map(i_q,:));
%     hold on;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust axes, legend properties etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ylim([0 .4]);


hLegend=legend('Location','SouthWest','Interpreter','latex',...
        'NumColumns',1);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')

hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$\omega_0/q$','interpreter','latex');
h_axis = gca;

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)


annotation_str = {modelname,sprintf('$N = (%d)^2$',sqrtN)};
dim=[.165 .17+.04*N_q .1 .2];
h_annotation=annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

figname=sprintf('%s/%s_omega_0_byq_vs_T',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end


%% 5.1 Setting it to omega_0 vs T

for i_q = 1:N_q
    q = q_vals_collect(i_q);
    h_fit_data(i_q).YData = h_fit_data(i_q).YData .* q;
    h_sim_data(i_q).YData = h_sim_data(i_q).YData .* q;
end
ylim([0 .17]);
set(hYLabel,'String','$\omega_0$');
set(hLegend,'Location','NorthWest');
% set(gca, 'Xscale', 'lin','Yscale', 'lin')
dim=[.43 .73 .1 .2];
set(h_annotation,'Position',dim,'FitBoxToText','on')

figname=sprintf('%s/%s_omega_0_vs_T',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

























%% 6 Plot omega_1^2 and gamma^2/4 vs T for a few q
figure
c_map=linspecer(N_q);
markersym_gamma = 's';
markersym_omega = '^';
for i_q = 1:N_q
    q = q_vals_collect(i_q);
    
    hold on;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot for omega_1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dispname = sprintf('$q = %.2f, \\omega_1$',q);
    dispname = sprintf('$q = %.2f$',q);
    data_cur = omega_1_vals(:,i_q).^2 ./ gamma_vals(:,i_q).^2 * 4;
    h_fit_data(i_q) = plot(T_vals,data_cur,...
        'LineStyle', '-.', ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_q,:));

    h_sim_data(i_q) = plot(T_vals,data_cur,'s',...
        'LineStyle', 'none', ...
        'Marker', markersym_omega, 'MarkerSize', 3.8, ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_q,:),...
        'MarkerFaceColor',c_map(i_q,:),'MarkerEdgeColor','k');

    h_leg_data(i_q) = plot(NaN,NaN,...
        'LineStyle', '-.', ...
        'LineWidth', 1.3, ...
        'Marker', markersym_omega, 'MarkerSize', 3.8, ...
        'DisplayName', dispname, ...
        'MarkerFaceColor',c_map(i_q,:),'MarkerEdgeColor','k', ...
        'Color',c_map(i_q,:));

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Plot for omega_1
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dispname = sprintf('$q = %.2f, \\gamma$',q);
%     data_cur = gamma_vals(:,i_q).^2 / 4;
%     h_fit_data(i_q) = plot(T_vals,data_cur,...
%         'LineStyle', '--', ...
%         'LineWidth', 2, ...
%         'HandleVisibility', 'off', ...
%         'Color',c_map(i_q,:));
% 
%     h_sim_data(i_q) = plot(T_vals,data_cur,'s',...
%         'LineStyle', 'none', ...
%         'Marker', markersym_gamma, 'MarkerSize', 3.8, ...
%         'LineWidth', 1, ...
%         'HandleVisibility', 'off', ...
%         'Color',c_map(i_q,:),...
%         'MarkerFaceColor',c_map(i_q,:),'MarkerEdgeColor','k');
% 
%     h_leg_data(i_q) = plot(NaN,NaN,...
%         'LineStyle', '--', ...
%         'LineWidth', 1.3, ...
%         'Marker', markersym_gamma, 'MarkerSize', 3.8, ...
%         'DisplayName', dispname, ...
%         'MarkerFaceColor',c_map(i_q,:),'MarkerEdgeColor','k', ...
%         'Color',c_map(i_q,:));
end

yline(1,'-','Color','k','HandleVisibility','off')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust axes, legend properties etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ylim([1e-2 Inf]);

hLegend=legend('Location','SouthWest','Interpreter','latex',...
        'NumColumns',1);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'YGrid', 'on', 'XGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', .5,'Xscale', 'lin','Yscale', 'log')

hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$4\omega_1^2 / \gamma^2$','interpreter','latex');
h_axis = gca;

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)


annotation_str = {modelname,sprintf('$N = (%d)^2$',sqrtN)};
dim=[.165 .17+.04*N_q .1 .2];
h_annotation=annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

figname=sprintf('%s/%s_omega_1_gamma_vs_T',basedir,curmodel);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end