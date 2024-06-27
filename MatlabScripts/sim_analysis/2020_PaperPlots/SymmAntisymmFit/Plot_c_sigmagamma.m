clear
cd /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/SymmAntisymmFit
run ../initialization_script
saveswitch=1;
basedir=sprintf('%s/SymmAntisymmFit',fig_base);

sqrtN_vals = [16 32 64 128];
sampfilename="samp_Dynamics";
    
    
% model="mxy"; "xy_s"
model="mxy";
if model=="mxy"
    curmodel="mxy";
    modeltag="MXY";
    sqrtN_vals = [16 32 64 128 256];
    L_vals=[9.25,18.5,37,74,148];
    T_vals = [.03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25];
    T_dirs = {"T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25"};
    dir_lintime="/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00";
    dir_reduced_dt="/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00";
    dir_longer_t="/data/scc/thobi/211201_LongerTime/mxy_3.00";
    simdir_helicity = "210715_LinearTimeSampling";
    modeldir="mxy_3.00";

    xlim_vec=[0 .25];
    ylim_sigma_vec=[0 3];
    ylim_omega=[0 .3];
    xlim_omega_inset=[0 .19];
    ylim_omega_inset=[.97 1.03];
    YTick_Upsilon_vec=(0:.01:2);

    T_BKT=.173;

elseif model=="xy_s"
    modeltag="XY";
    model="xy_s";
    curmodel="xy";
    sqrtN_vals = [16 32 64 128];
    L_vals = sqrtN_vals;
    T_dirs = {"T_.10", "T_.20", "T_.30", "T_.40", "T_.50", "T_.60", "T_.70", "T_.80", "T_.85", "T_.91", "T_.95", "T_1.00"};
    T_vals = [.10, .20, .30, .40, .50, .60, .70, .80, .85, .91, .95, 1.00];
    dir_reduced_dt="/data/scc/thobi/220201_ReducedSmapleStepDeltat/xy_s";
    simdir_helicity = "220201_ReducedSmapleStepDeltat";
    modeldir="xy_s";

    xlim_vec=[.5 1];
    ylim_sigma_vec=[0 3];
    ylim_omega=[0 1];
    xlim_omega_inset=[.5 .9];
    ylim_omega_inset=[.89 1.11];
    YTick_Upsilon_vec=.9:.05:1.1;

end
rho_vals = sqrtN_vals.^2 ./ L_vals .^2;
savefile=sprintf('%s/%s_dataset',basedir,model);
%load(sprintf('%s/%s_dataset',basedir,model));

fitfunc_symm="exp(-gamma*x/2)*(cos(omega_1*x) + gamma/2/omega_1*sin(omega_1*x))";
fitfunc_asymm="exp(-gamma*x/2)*cos(omega_1*x)";
f_symm=@(t,om,ga) exp(-ga*t/2).*(cos(om*t)+ga/2/om*sin(om*t));




%% Fixed N, sigma error bars, only one fit
if model == "mxy"
%     T_select=[.11, .14, .16, .17, .175, .18, .185, .20];
    T_select=[.03, .07, .09, .11, .14, .16, .17, .175, .18, .185, .20, .22, .24];
    sqrtN_select=[16,32,64,128];
    sqrtN_select=[128,256];
    L_select=L_vals(4:5);
    
elseif model == "xy_s"
    T_select=[.50 .60 .70 .80 .85 .91 .95 1.00];
    sqrtN_select=[32,64,128];
    
end
L_select=L_vals(ismember(sqrtN_vals,sqrtN_select));
n_N=numel(sqrtN_select);
n_T=numel(T_select);
fittype="lin";
qmode="qint";

q_cell=cell(n_N,n_T);
fitob_ga_cell=cell(n_N,n_T);
file_cell=cell(n_N,n_T);
data_cell=cell(n_N,n_T);
ga_vals_cell=cell(n_N,n_T);
sigma_vals=zeros(n_N,n_T);
sigma_err_vals=zeros(n_N,n_T);
sigma_cimin_vals=zeros(n_N,n_T);
sigma_cimax_vals=zeros(n_N,n_T);

fitob_om_cell=cell(n_N,n_T);
om_1_vals_cell=cell(n_N,n_T);
c_vals=zeros(n_N,n_T);
c_err_vals=zeros(n_N,n_T);
c_cimin_vals=zeros(n_N,n_T);
c_cimax_vals=zeros(n_N,n_T);

H_x_vals=zeros(n_N,n_T);
H_y_vals=zeros(n_N,n_T);
I_x_2_vals=zeros(n_N,n_T);
I_y_2_vals=zeros(n_N,n_T);
Helicity_vals=zeros(n_N,n_T);
for i_N=1:n_N
    fprintf("N = (%d)^2\n",sqrtN_select(i_N));
    for i=1:numel(T_select)
        T=T_select(i);
        i_T=find(T_vals==T);
        sqrtN=sqrtN_select(i_N);
        T=T_vals(i_T);
        fprintf('T = %.3f ',T);

        
        curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_reduced_dt,sqrtN_select(i_N),T_dirs{i_T},sampfilename);
        q_off=1;
        if ~isfile(curfile)
            curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_longer_t,sqrtN_select(i_N),T_dirs{i_T},sampfilename);
            q_off=1;
            if ~isfile(curfile)
                curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_lintime,sqrtN_select(i_N),T_dirs{i_T},sampfilename);
                q_off=0;
            end
        end
        file_cell{i_N,i}=curfile;
        data_cell{i_N,i}=load(curfile,"averaging_times","gmperpmperp","qbin");
        t=data_cell{i_N,i}.averaging_times;
        qbin=data_cell{i_N,i}.qbin;
        gmperpmperp=data_cell{i_N,i}.gmperpmperp;
        
        if qmode == "qint"
            q_integers=find((qbin(1:end-q_off)/qbin(1)>round(qbin(1:end-q_off)/qbin(1))-1e-3) ...
                .* (qbin(1:end-q_off)/qbin(1)<round(qbin(1:end-q_off)/qbin(1))+1e-3));
        elseif qmode == "qfull"
            q_integers=1:numel(qbin)-q_off;
        end
        q_integers=find(gmperpmperp(q_integers)~=0);
        
        n_q = numel(qbin);
        n_q_int = numel(q_integers);
        ga_vals=zeros(1,n_q_int);
        om_vals=zeros(1,n_q_int);
        
        q_vals=qbin(q_integers);
        for i_q_int = 1:numel(q_vals)
            i_q = q_integers(i_q_int);
            fprintf('q = %.3f, ',q_vals(i_q_int));
            cf=real(gmperpmperp(i_q:numel(qbin):end))/real(gmperpmperp(i_q));
            fitob_sym=fit(t(:),cf(:),fitfunc_symm,'StartPoint',[1e-2, 1e-2], 'Lower', [0, 0]);
            ga_vals(i_q_int)=fitob_sym.gamma;
            om_vals(i_q_int)=abs(fitob_sym.omega_1);
        end
        q_cell{i_N,i}=q_vals;
        ga_vals_cell{i_N,i}=ga_vals;
        om_1_vals_cell{i_N,i}=om_vals;
        if fittype == "lin"
            fitfunc_pow="a*x^b";
            [fitob_ga,gof,out]=fit(q_vals(:),ga_vals(:),fitfunc_pow, ...
                'StartPoint',[ga_vals(1)/q_vals(1)^2,2], ...
                'Lower',[0,0]);
        elseif fittype == "log"
            fitfunc_pow="b*log(x)+a";
            [fitob_ga,gof,out]=fit(q_vals(:),log(ga_vals(:)),fitfunc_pow, ...
                'StartPoint',[log(ga_vals(1)/q_vals(1)^2),2], ...
                'Lower',[-Inf,0]);
            
        end
        
        fitob_ga_cell{i_N,i}=fitob_ga;
        sigma_vals(i_N,i)=fitob_ga.b;

        fitfunc_c="c*x";
%         om_select=find(ga_vals<om_vals);
        om_select=1:4;
        if isempty(om_select)
            c_vals(i_N,i)=om_vals(1) / q_vals(1);
        else
            [fitob_om,gof,out]=fit(q_vals(om_select)',om_vals(om_select)',fitfunc_c, ...
                'StartPoint',[1], ...
                'Lower',[0,0]);
            fitob_om_cell{i_N,i}=fitob_om;
            c_vals(i_N,i)=fitob_om.c;
        end
        
        ci=confint(fitob_ga);
        sigma_cimin_vals(i_N,i)=ci(1,2);
        sigma_cimax_vals(i_N,i)=ci(2,2);
        sigma_err_vals(i_N,i)=(ci(2,2)-ci(1,2))/3.92;
    
        ci=confint(fitob_om);
        c_cimin_vals(i_N,i)=ci(1);
        c_cimax_vals(i_N,i)=ci(2);
        c_err_vals(i_N,i)=(ci(2)-ci(1))/3.92;
        
        
        filename_Helicity=sprintf('/data/scc/thobi/%s/%s/sqrtN_%d/%s/samp_Dynamics_helicity.mat',...
                simdir_helicity,modeldir,sqrtN,T_dirs{i_T});
        S_helicity=load(filename_Helicity,...
                "H_x","H_y","I_x_2","I_y_2");
        H_x_vals(i_N,i)=S_helicity.H_x;
        H_y_vals(i_N,i)=S_helicity.H_y;
        I_x_2_vals(i_N,i)=S_helicity.I_x_2; 
        I_y_2_vals(i_N,i)=S_helicity.I_y_2;
        Helicity_vals(i_N,i)=1/2/L_select(i_N)^2 * ...
            ( abs(H_x_vals(i_N,i) + H_y_vals(i_N,i)) - ...
            1/T_select(i)*(I_x_2_vals(i_N,i) + I_y_2_vals(i_N,i)));
        fprintf("\n");


    end
    
end
save(savefile);
%%
load(savefile);
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
set(gca,'FontSize',fontsize_axis);
c_map=colormap_sqrtN();
for i_N=1:n_N
    sqrtN=sqrtN_select(i_N);
    i_N_full=find(sqrtN_vals == sqrtN);
    dispname="$N =(" + sqrtN + ")^2$";
%     plot(T_select,sigma_cimin_vals(i_N,:),'xr','DisplayName','HandleVisbility','off'); hold on;
%     plot(T_select,sigma_cimax_vals(i_N,:),'xg','DisplayName','ci max'); 
    errorbar(T_select,sigma_vals(i_N,:),sigma_err_vals(i_N,:), ...
        '-','DisplayName',dispname, ...
        'Color',c_map(i_N_full,:)); hold on;
end
hlegend=legend('Location','best','Interpreter','latex',...
    'FontName','cmr14','fontsize',fontsize_labels);
hxlabel=xlabel('$T$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
hylabel=ylabel('$\sigma_{\gamma}$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
htitle=title(sprintf('$N=(%d)^2$',sqrtN),'Interpreter','latex',...
    'FontName','cmr14','fontsize',fontsize_titles);

subtitle_txt=sprintf('Fittype "%s", q selection "%s"',fittype,qmode);
hsubtitle=subtitle(subtitle_txt,'Interpreter','latex',...
    'FontName','cmr14','fontsize',fontsize_labels);


figname=sprintf('%s/figs/%s_ga_fullT_err_%sfit_%s_sqrtN_%d',basedir,curmodel,fittype,qmode,sqrtN);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');


%% 0.75 Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall settings, annotatiosn etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
load(savefile);
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm 1.8*columnwidth_cm]);
linestyle='-';
linewidth=1;
markersize=5;

% colors={'r', 'b'};
markers={'s', 'd', 'o', '^','v'};


y_border = .008;
x_border = .008;
labelbox_size = [.05 .03];


c_map=colormap_sqrtN();








    
%% 2 Subplot bottom: sigma_gamma fit vs FFT

% subplot(2,1,2,'replace');
% inset_pos=[pos(1)+.005,...
%     pos(2)+.4*pos(4),...
%     .6*pos(3),...
%     .56*pos(4)];
% ax_inset = axes('OuterPosition',inset_pos);
% ax_bot=axes('InnerPosition',[.2 .1 .75 .4]);
% pos_plotbottom=get(gca,'Position');
ax_bot=subplot(2,1,2,'replace');
set(ax_bot,'InnerPosition',[.2 .1 .75 .4]);
pos_plottop=get(gca,'Position');
hold on;



for i_N = 1:n_N
    sqrtN=sqrtN_select(i_N);
    i_N_full=find(sqrtN_vals == sqrtN);
    dispname=sprintf('$N = (%d)^2$',sqrtN);
    y=sigma_vals(i_N,:);
    y_err=sigma_err_vals(i_N,:);
    h_plot{i_N} = errorbar(T_select,y,y_err,...
        'DisplayName',dispname,...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'Marker', markers{i_N},...
        'MarkerSize', markersize,...
        'HandleVisibility', 'on', ...
        'Color',c_map(i_N_full,:));
end
xlim(xlim_vec);
ylim(ylim_sigma_vec);


hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$\sigma_{\gamma}$','interpreter','latex');
hLegend = legend('Location', 'SouthWest','interpreter','latex',...
    'NumColumns',1);
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','lin',...
    'XAxisLocation','bottom',...
    'LineWidth', .5)
% ylim([1e1 Inf]);

xlim_vec=get(gca,"xlim");
ylim_vec=get(gca,"ylim");

xline(T_BKT,'LineStyle','--','HandleVisibility','off');
text(.97*T_BKT,.08*ylim_vec(2),'$T_{\textrm{BKT}}$',...
        'VerticalAlignment','middle','HorizontalAlignment','right',...
        'interpreter','latex',...
        'Color','k','FontSize', fontsize_annotation);


% pos = get(gca,'Position');
% dim=[pos(1)+pos(3)-x_border-labelbox_size(1), ...
%     pos(2)+y_border, ...
%     labelbox_size(1),...
%     labelbox_size(2)];
% textpos=[.94 .06];
% text('Position',textpos, 'String', '(b)', ...
%     'interpreter','latex',...
%     'Units','normalized',...
%     'Color','black','FontSize', fontsize_subfiglabels);
h_text=add_subfig_label(gca,"$(b)$","nw","lin","lin",fontsize_subfiglabels);






%% 1 Subplot top: omega
ax_top=subplot(2,1,1,'replace');
set(ax_top,'InnerPosition',[.2 .5 .75 .4]);
pos_plottop=get(gca,'Position');
hold on;


c_map=colormap_sqrtN();
for i_N = 1:n_N
    sqrtN=sqrtN_select(i_N);
    i_N_full=find(sqrtN_vals == sqrtN);
    dispname=sprintf('$N = (%d)^2$',sqrtN);
%     y=c_sw_fit(i_N,:);
    y=c_vals(i_N,:);
    h_plot{i_N} = plot(T_select,y,...
        'DisplayName',dispname,...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'Marker', markers{i_N},...
        'MarkerSize', markersize,...
        'HandleVisibility', 'on', ...
        'Color',c_map(i_N_full,:));
end
xlim(xlim_vec);
ylim(ylim_omega);


hXLabel = xlabel('$T$','interpreter','latex');
% hYLabel = ylabel('$S_{m\perp}(q,\omega) / \chi_{m\perp}(q)$','interpreter','latex');
% hYLabel = ylabel('$\gamma,S^{\textrm{max}}$','interpreter','latex');
hYLabel = ylabel('$c$','interpreter','latex');
% hLegend = legend('Location', 'SouthWest','interpreter','latex',...
%     'NumColumns',1);
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
% set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','lin',...
    'XAxisLocation','top',...
    'LineWidth', .5)
% ylim([1e1 Inf]);


% pos = get(gca,'Position');
% dim=[pos(1)+pos(3)-x_border-labelbox_size(1), ...
%     pos(2)+y_border, ...
%     labelbox_size(1),...
%     labelbox_size(2)];
% textpos=[.94 .06];
% text('Position',textpos, 'String', '(a)', ...
%     'interpreter','latex',...
%     'Units','normalized',...
%     'Color','black','FontSize', fontsize_subfiglabels);
h_text=add_subfig_label(gca,"$(a)$","nw","lin","lin",fontsize_subfiglabels);


annotation_str = {modeltag};
pos = get(gca,'InnerPosition');
dim=[pos(1)+pos(3)-.14, ...
    pos(2)+pos(4)-.075, ...
    .125,.09];
h_annot=annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

xlim_vec=get(gca,"xlim");
ylim_vec=get(gca,"ylim");

xline(T_BKT,'LineStyle','--','HandleVisibility','off');
text(.97*T_BKT,.9*ylim_vec(2),'$T_{\textrm{BKT}}$',...
        'VerticalAlignment','middle','HorizontalAlignment','right',...
        'interpreter','latex',...
        'Color','k','FontSize', fontsize_annotation);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=get(gca,"InnerPosition");
axes("Position",[pos(1)+.08, pos(2)+.03, pos(3)/2, pos(4)/2])

for i_N = 1:n_N
    sqrtN=sqrtN_select(i_N);
    i_N_full=find(sqrtN_vals == sqrtN);
    dispname=sprintf('$N = (%d)^2$',sqrtN);
    h_plot{i_N} = plot(T_select,rho_vals(i_N)*c_vals(i_N,:).^2./Helicity_vals(i_N,:),...
        'DisplayName',dispname,...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'Marker', 'o',...
        'MarkerSize', markersize,...
        'HandleVisibility', 'on', ...
        'Color',c_map(i_N_full,:));
    hold on;
end
xlim(xlim_omega_inset);
ylim(ylim_omega_inset);
text(xlim_omega_inset(1)+0.01,ylim_omega_inset(2)-.002,'$\frac{c^2}{\Upsilon / \rho I}$',...
        'VerticalAlignment','top','HorizontalAlignment','left',...
        'interpreter','latex',...
        'Color','k','FontSize', fontsize_annotation);
text(xlim_omega_inset(2)-0.01,ylim_omega_inset(1)+.002,'$T$',...
        'VerticalAlignment','bottom','HorizontalAlignment','right',...
        'interpreter','latex',...
        'Color','k','FontSize', fontsize_annotation);    
set(gca, 'FontName', 'cmr12','FontSize', fontsize_axis,...
    'YTick',YTick_Upsilon_vec);


















%% Saving
figname=sprintf('%s/%s_omegagamma_SpinwaveGammaexp',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end



%% 3 omega no subplot
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
pos_plottop=get(gca,'Position');
hold on;


c_map=colormap_sqrtN();
for i_N = 1:n_N
    sqrtN=sqrtN_select(i_N);
    i_N_full=find(sqrtN_vals == sqrtN);
    dispname=sprintf('$N = (%d)^2$',sqrtN);
%     y=c_sw_fit(i_N,:);
    y=c_vals(i_N,:);
    h_plot{i_N} = plot(T_select,y,...
        'DisplayName',dispname,...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'Marker', markers{i_N},...
        'MarkerSize', markersize,...
        'HandleVisibility', 'on', ...
        'Color',c_map(i_N_full,:));
end
xlim(xlim_vec);
ylim(ylim_omega);


hXLabel = xlabel('$T$','interpreter','latex');
% hYLabel = ylabel('$S_{m\perp}(q,\omega) / \chi_{m\perp}(q)$','interpreter','latex');
% hYLabel = ylabel('$\gamma,S^{\textrm{max}}$','interpreter','latex');
hYLabel = ylabel('$c$','interpreter','latex');
% hLegend = legend('Location', 'SouthWest','interpreter','latex',...
%     'NumColumns',1);
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
% set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','lin',...
    'XAxisLocation','bottom',...
    'LineWidth', .5)


annotation_str = {modeltag};
pos = get(gca,'InnerPosition');
dim=[pos(1)+pos(3)-.14, ...
    pos(2)+pos(4)-.075, ...
    .125,.09];
h_annot=annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=get(gca,"InnerPosition");
axes("Position",[pos(1)+.08, pos(2)+.05, pos(3)/2, pos(4)/2])

for i_N = 1:n_N
    sqrtN=sqrtN_select(i_N);
    i_N_full=find(sqrtN_vals == sqrtN);
    dispname=sprintf('$N = (%d)^2$',sqrtN);
    h_plot{i_N} = plot(T_select,c_vals(i_N,:).^2./Helicity_vals(i_N,:),...
        'DisplayName',dispname,...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'Marker', 'o',...
        'MarkerSize', markersize,...
        'HandleVisibility', 'on', ...
        'Color',c_map(i_N_full,:));
    hold on;
end
xlim(xlim_omega_inset);
ylim(ylim_omega_inset);
text(xlim_omega_inset(1)+0.02,ylim_omega_inset(2)-.002,'$c^2 / \Upsilon$',...
        'VerticalAlignment','top','HorizontalAlignment','left',...
        'interpreter','latex',...
        'Color','k','FontSize', fontsize_annotation);
text(xlim_omega_inset(2)-0.01,ylim_omega_inset(1)+.002,'$T$',...
        'VerticalAlignment','bottom','HorizontalAlignment','right',...
        'interpreter','latex',...
        'Color','k','FontSize', fontsize_annotation);    
set(gca, 'FontName', 'cmr12','FontSize', fontsize_axis,...
    'YTick',YTick_Upsilon_vec);




figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/Magnetization/SpinACF_Fit/%s_c_plot',model);
% figname=sprintf('%s/%s_omegagamma_SpinwaveGammaexp',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end





















%% 4 sigma_gamma fit vs FFT, no subplot
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

pos_plottop=get(gca,'Position');
hold on;


c_map=colormap_sqrtN();
for i_N = 1:n_N
    sqrtN=sqrtN_select(i_N);
    i_N_full=find(sqrtN_vals == sqrtN);
    dispname=sprintf('$N = (%d)^2$',sqrtN);
    y=sigma_vals(i_N,:);
    y_err=sigma_err_vals(i_N,:);
    h_plot{i_N} = errorbar(T_select,y,y_err,...
        'DisplayName',dispname,...
        'LineStyle', linestyle, ...
        'LineWidth', linewidth, ...
        'Marker', markers{i_N},...
        'MarkerSize', markersize,...
        'HandleVisibility', 'on', ...
        'Color',c_map(i_N_full,:));
end
xlim(xlim_vec);
ylim(ylim_sigma_vec);


hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$\sigma_{\gamma}$','interpreter','latex');
hLegend = legend('Location', 'SouthWest','interpreter','latex',...
    'NumColumns',1);
h_axis = gca;
% Font
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'YScale','lin', 'XScale','lin',...
    'XAxisLocation','bottom',...
    'LineWidth', .5)



figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/Magnetization/SpinACF_Fit/%s_sigmagamma',model);
% figname=sprintf('%s/%s_omegagamma_SpinwaveGammaexp',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    %                 print(sprintf('%s.eps',figname),'-depsc');
end
