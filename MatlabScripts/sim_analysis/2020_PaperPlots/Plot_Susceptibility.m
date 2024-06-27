%% 0 Initialization
run initialization_script;
basedir=sprintf('%s/plots/static',fig_base);

saveswitch=1;
i_model = 1;
if (i_model == 1)
    curmodel="mxy";
    curtitle="MXY model";
    dir='/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00';
    sqrtN_vals=[16 32 64 128 256];
    T_dirs = {"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40"};
    T_vals=[.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40];
    T_dirs = {"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33"};
    T_vals=[.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33];

    collectfilename='samp_Dynamics';
    collectfilename='samp_mag_errors';
%     data=load(mxydata_name);
    runmax=500;

    L_vals=[9.25,18.5,37,74,148];
elseif (i_model == 2)
    curmodel="xy";
    curtitle="SXY model";
    collectfilename='samp_Dynamics';
    
    data=load(xydata_name);
    runmax=250;

    L_vals=[16,32,64,128,256];
    T_max = 3.5;
elseif (i_model == 3)
    curmodel="fmxy";
    curtitle="FMXY model"; 
    collectfilename='samp_Dynamics_collect';
    
    data=load(fmxydata_name);
    runmax=250;
    L_vals=[9.25,18.5,37,74,148];
    T_max = .4;
elseif (i_model == 4)
    curmodel="xy_s";
    curtitle="SXY model";
    runmax=125;
    collectfilename='samp_eq_collect';
    
    data=matfile(xysdata_name);
    FS_Mag_fitdata=matfile(xysfit_FSMag_name);

    L_vals=[16 32 64 128 256];
    T_max=2;
end
N_N = numel(sqrtN_vals);
N_T = numel(T_dirs);
N_vals=sqrtN_vals.^2;
savefile=sprintf('data_%s_Susceptibility',curmodel);
T_scaling = "off";

%% 1 Data Assembly
absM_mean=zeros(N_N,N_T);
absM_err=zeros(N_N,N_T);
M_2_mean=zeros(N_N,N_T);
M_2_err=zeros(N_N,N_T);
M_4_mean=zeros(N_N,N_T);
M_4_err=zeros(N_N,N_T);
sigma_vals=absM_mean;
sigma_err_vals=absM_mean;
chi_mean=zeros(N_N,N_T);
chi_err=zeros(N_N,N_T);
sigmamax_mean=zeros(N_N,1);
sigmamax_err=zeros(N_N,1);
sigmamax_i_T=zeros(N_N,1);
sigmamax_T=zeros(N_N,1);
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    N = sqrtN^2;
    fprintf('sqrtN = %d',sqrtN);
    for i_T = 1:N_T
        T = T_vals(i_T);
        fprintf(' %.3g',T);
        curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir,sqrtN,T_dirs{i_T},collectfilename);
        load(curfile,'m_mean','mvar_mean','mvar_err');
        absM_mean(i_N,i_T)=m_mean;
        sigma_vals(i_N,i_T)=mvar_mean;
        sigma_err_vals(i_N,i_T)=mvar_err;
%         load(curfile,'absM_av','M_2_av','M_4_av')
%         absM_mean(i_N,i_T)=absM_av;
%         M_2_mean(i_N,i_T)=M_2_av;
%         M_4_mean(i_N,i_T)=M_4_av;
    end
%     sigma_vals(i_N,:) = M_2_mean(i_N,:) - absM_mean(i_N,:).^2;
    chi_mean(i_N,:) = N * sigma_vals(i_N,:);
    if T_scaling == "on"
        chi_mean(i_N,i_T)=chi_mean(i_N,:)./T_vals;
    end

    [sigmamax,i_T_max]=max(sigma_vals(i_N,:));
    
    T = T_vals(i_T_max);
    sigmamax_T(i_N) = T;
    sigmamax_i_T(i_N) = i_T_max;
    fprintf('\n   Max at T=%.3g\n',T);
%     curfile=sprintf('%s/sqrtN_%d/%s/%s_collect.mat',dir,sqrtN,T_dirs{i_T_max},collectfilename);
%     load(curfile,'M_2_av_collect','absM_av_collect');
%     chicur=N*(M_2_av_collect - absM_mean(i_N,i_T_max).^2);
    sigmamax_mean(i_N)= sigmamax;%mean(chicur);
    sigmamax_err(i_N)= sigma_err_vals(i_N,i_T_max);

%     varerr_vec=sqrt(1/runmax * (M_2_mean(i_N,i_T_max) - (runmax - 3) / (runmax - 1) *var_vec.^2));
%     chierr{i_N} = N_vals(i_N) * varerr_vec ./ T_vals;
%     chierr_vec(i_N) = N_vals(i_N) * sqrt(1/runmax * (M_2_mean(i_N,i_T_max) ...
%         - (runmax - 3) / (runmax - 1) *(M_2_mean(i_N,i_T_max) - absM_mean(i_N,i_T_max).^2).^2));
    
%     if T_scaling == "on"
%         chierr_vec(i_N) = chierr_vec(i_N) / T;
%     end

%     var_vec = M_2_av(i_N,:) - absM_av(i_N,:).^2;
%     chiabsm_vec = N_vals(i_N) * var_vec ./ T_vals;
%     chiabsm{i_N} = chiabsm_vec;
%     
%     varerr_vec=sqrt(1/runmax * (M_2_av(i_N,:) - (runmax - 3) / (runmax - 1) *var_vec.^2));
%     % Apparently pg 483, Rao, 1973
%     chierr{i_N} = N_vals(i_N) * varerr_vec ./ T_vals;
end
% sigma_var_mat = M_4_mean - M_2_mean.^2;
% sigma_err_vals = 0 *sqrtN_vals;
% for i_N = 1:numel(sqrtN_vals)
%     sigma_err_vals(i_N) = 1/sqrt(runmax)*sqrt(sigma_var_mat(i_N,sigmamax_i_T(i_N)));
% end
%% Fit
eta_fitobs_absM=cell(1,N_T);
eta_vals_absM=0*T_vals;
eta_fitobs_chi=cell(1,N_T);
eta_vals_chi=0*T_vals;
for i_T = 1:N_T
    T = T_vals(i_T);
    absM_vec=absM_mean(:,i_T);
    eta_fitobs_absM{i_T}=fit_eta_Magnetization_FS(absM_vec,sqrtN_vals);
    eta_vals_absM(i_T)=eta_fitobs_absM{i_T}.eta;

    chi_vec=T*chi_mean(:,i_T);
    eta_fitobs_chi{i_T}=fit_eta_Susceptibility_FS(chi_vec,sqrtN_vals);
    eta_vals_chi(i_T)=eta_fitobs_chi{i_T}.eta;
end
save(savefile);



%% 2 The plot
load(savefile);
figure
hold on
c_map = colormap_sqrtN();
N_vals=sqrtN_vals.^2;
for i_N = 1 : N_N
    sqrtN=sqrtN_vals(i_N);
    N=sqrtN^2;
    chiabsm_vec=chi_mean(i_N,:);

%     exponents=2*(1-.5*eta_vals_absM);
    exponents=0;
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    hchiabsM(i_N) = plot(T_vals,chiabsm_vec./sqrtN_vals(i_N).^exponents,...
        'LineStyle', '-', ...
        'LineWidth', .8,...
        'Marker', 's', 'MarkerSize', 4, ...
        'DisplayName', dispname, ...
        'Color',c_map(i_N,:));
%     hchiabsM_line(i_N) = line(T_vals,chiabsm_vec./sqrtN_vals(i_N).^exponents);
%     set(hchiabsM_line(i_N), ...
%         'LineStyle', '-', 'LineWidth', 1, ...
%         'HandleVisibility','off',...
%         'Color', c_map(i_N,:))
%         hchiabsM_dots(i_N) = line(T_vals,chiabsm_vec);
%     dispname = sprintf('$N = (%d)^2$',sqrtN);
%     hchiabsM_dots(i_N) = line(T_vals,chiabsm_vec./sqrtN_vals(i_N).^exponents);
%     set(hchiabsM_dots(i_N), ...
%         'LineStyle', 'none', ...
%         'Marker', 'o', 'MarkerSize', 5, ...
%         'HandleVisibility', 'off', ...
%         'Color',c_map(i_N,:),...
%         'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:));
%     
%     hlabel(i_N)=plot(NaN,NaN,...
%         'LineStyle', '-', ...
%         'LineWidth', 1,...
%         'Marker', 'o', 'MarkerSize', 5, ...
%         'DisplayName', dispname, ...
%         'Color',c_map(i_N,:),...
%         'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:));
    
end
xlim([0, max(T_vals)]);
ylim([0 Inf]);
% xlim([0 .2]);
% ylim([0 .16]);


% hLegend = legend('Location', 'NorthWest','interpreter','latex');

hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$\chi_{m}$','interpreter','latex');
h_axis = gca;
% hTitle = title(curtitle);

% 3a Font
set(gca, 'FontSize', fontsize_axis, 'FontName', 'cmr12');
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
% set_fonts_default


% 3a Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', ...
    'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', .5)
ax_full = gca;

% 3a Inset
% ax_inset = axes('Position',[.62 .62 .26 .26]);
pos=get(ax_full,'InnerPosition');
ax_inset_pos = [pos(1)+.08,...
    pos(2)+.6*pos(4),...
    .35*pos(3),...
    .35*pos(4)];

% ax_inset_pos = [pos(1)+.6*pos(3),...
%     pos(2)+.5*pos(4),...
%     .35*pos(3),...
%     .45*pos(4)];
ax_inset = axes('Position',ax_inset_pos);

for i_N = 1:N_N
    sqrtN=sqrtN_vals(i_N);
    N=sqrtN^2;
    chiabsm_vec=chi_mean(i_N,:);
    exponents=2*(1-.5*eta_vals_absM);
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    hchiscaled(i_N) = plot(T_vals,chiabsm_vec./sqrtN_vals(i_N).^exponents,...
        'LineStyle', '-', ...
        'LineWidth', .5,...
        'Marker', 's', 'MarkerSize', 3, ...
        'DisplayName', dispname, ...
        'Color',c_map(i_N,:));
    hold on;
end

xlim([0 .2]);
ylim([0 .01]);
set(ax_inset, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'off', 'XGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'XTick', 0:.05:2, 'YTick', 0:.002:1, ...
    'LineWidth', .5, 'Xscale', 'lin','Yscale','lin')
set(gca, 'FontSize', fontsize_axis, 'FontName', 'cmr12');

% 3a Text
NW = [min(xlim) min(ylim)+.98*diff(ylim)];
SE = [max(xlim) min(ylim)];
text(SE(1),SE(2),'$T$',...
    'VerticalAlignment','bottom','HorizontalAlignment','right',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_axis);
text(NW(1),NW(2),'\ $\chi_m/ N^{1-\eta/2}$',...
    'VerticalAlignment','top','HorizontalAlignment','left',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_axis); %fontsize_ax_labels

% Second inset

pos=get(ax_full,'InnerPosition');
ax_inset_r_pos = [pos(1)+.7*pos(3),...
    pos(2)+.6*pos(4),...
    .28*pos(3),...
    .35*pos(4)];
ax_inset_right = axes('Position',ax_inset_r_pos);

% 
% xlim_vec=[.5*N_vals(1),2*N_vals(end)];
% ylim_vec=[0 .03];
% xlim(xlim_vec)
% ylim(ylim_vec)

% chimax_err chimax_mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hchimax_line = errorbar(N_vals,chimax_mean(:)./N_vals(:),chimax_err(:)./N_vals(:));
f_POW = fittype('a*(x)^(-eta)');
init_exponent = .25;
init_factor = 1;
c0=[init_factor,init_exponent];
Lower_lim=[0, 0];
Upper_lim=[Inf, Inf];

%     c = fit(sqrtN_vals(1:end-1)', chimax_mean(1:end-1)', f_POW,'StartPoint',c0,...
%         'Lower',Lower_lim, 'Upper', Upper_lim);
c = fit(sqrtN_vals(1:end)', sigmamax_mean(1:end), f_POW,'StartPoint',c0,...
    'Lower',Lower_lim, 'Upper', Upper_lim,...
    'Algorithm','Trust-Region', 'Robust', 'LAR');

c_log = fit(log(sqrtN_vals(1:end).^2)', log(sigmamax_mean(1:end)), "-a*x + b",'StartPoint',c0,...
    'Algorithm','Trust-Region', 'Robust', 'LAR');
eta_log = 2*c_log.a;

f_POW_fixedeta = fittype('a*(x)^(-.25)');
init_factor = 1;
c0=[init_factor];
Lower_lim=[0];
Upper_lim=[Inf];

c_fixedeta = fit(sqrtN_vals(:), sigmamax_mean(:), f_POW_fixedeta,'StartPoint',c0,...
    'Lower',Lower_lim, 'Upper', Upper_lim);



%     hchimax_fit = line(NN,c.a * sqrt(NN).^(- c.eta));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on;

NN=logspace(1,6);
% hchimax_fit = plot(NN,c_fixedeta.a * sqrt(NN).^(- .25),...
% hchimax_fit = plot(NN,c.a * sqrt(NN).^(- c.eta),...
%     'LineStyle','-',...
%     'LineWidth',1,...
%     'Color','b');
hchimax_fit = plot(NN,c_fixedeta.a * sqrt(NN).^(- .25),...
    'LineStyle','-',...
    'LineWidth',1,...
    'Color','b');
% hchimax_line = errorbar(N_vals,sigmamax_mean'./N_vals,sqrt(chierr_vec)./N_vals,... Factor 1/2 because above and below
hchimax_line = errorbar(N_vals,sigmamax_mean,sigmamax_err,... 
    'LineWidth',1,...
    'LineStyle','none',...
    'Marker','x','MarkerSize',4,...
    'CapSize', 3,...
    'Color','r');
% hchimax_line = plot(N_vals,sigmamax_mean,sigma_err_vals,... Factor 1/2 because above and below
%     'LineWidth',1,...
%     'LineStyle','none',...
%     'Color','r');
% set(hchimax_line, ...
%     'LineStyle', 'none', 'LineWidth', 1, ...
%     'Marker', 'none', 'MarkerSize', 5, ...
%     'Color',c_map(1,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(1,:))

xlabelstr=sprintfc("(%d)^2",sqrtN_vals);
set(ax_inset_right, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'off', 'XGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'XTick', 16^2*4.^(0:4), 'YTick', 0:1e-2:5e-2, ...
    'XTickLabel',xlabelstr,...
    'LineWidth', .5, 'Xscale', 'log','Yscale','lin')
% % xticklabels()
% 
% % hXLabel = xlabel('$N$','interpreter','latex');
% % hYLabel = ylabel('$\chi_{\textrm{max}} / N$','interpreter','latex');
% 
% % set_fonts_default
set(gca, 'FontSize', fontsize_axis, 'FontName', 'cmr12');

xlim_vec=[.5*N_vals(1),4*N_vals(end)];
ylim_vec=[0 .031];
xlim(xlim_vec)
ylim(ylim_vec)

% 3a Text
% NW = [min(xlim) max(ylim)]+[log(diff(xlim)) -log(diff(ylim))]*0.02;
% SE = [max(xlim) min(ylim)]+[-diff(xlim) .5*diff(ylim)]*0.02;
pos=get(gca,'Position');
% NW=[pos(1),pos(2)+pos(4),.1,.1];
% SE=[pos(1)+pos(3),pos(2),.1,.1];
% h_ylabel=annotation("textbox",NW,'String',"$\chi_{\textrm{max}} / N$",...
%     "Units","normalized","FitBoxToText","on",...
%     "Interpreter","latex",...
%     "EdgeColor","none",...
%     "HorizontalAlignment","left",...
%     "VerticalAlignment","top");
% pos=get(h_ylabel,"Position");
% set(h_ylabel,"Position",[pos(1), pos(2)-pos(4), pos(3), pos(4)]);
NW = [min(xlim) min(ylim)+.98*diff(ylim)];
SE = [max(xlim) min(ylim)];
SW = [min(xlim) min(ylim)];
text(SE(1),SE(2),'$N$',...
    'VerticalAlignment','bottom','HorizontalAlignment','right',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_axis);
text(NW(1),NW(2),'\ $\sigma_{\textrm{max}}^2$',...
    'VerticalAlignment','top','HorizontalAlignment','left',...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_axis);

text(SW(1),SW(2),sprintf('\\ $\\eta = %.2f$', .25),... c.eta),...
    'VerticalAlignment','bottom','HorizontalAlignment','left',...
    'interpreter','latex',...
    'Color','b','FontSize', fontsize_ax_labels);

% 
%  annotation('textbox','String',sprintf('$\\eta = %.2f$',.25),...
%      'Position',ax_inset.Position + [.07 .155 .0 .0],...
%      'LineWidth',.5,...
%      'HorizontalAlignment','left','Vert','bottom',...
%      'interpreter','latex','FitBoxToText','on',...
%      'BackgroundColor','White',...
%      'FontSize', fontsize_annotation, 'FontName', 'cmr12')

 set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

figname=sprintf('%s/%s_AbsSusceptibility',basedir,curmodel);
fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end





















