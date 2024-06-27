%% 0 Initialization
cd /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/NF_FFT
clear
run ../initialization_script;
saveswitch=1;
basedir="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/NF_FFT";


curmodel="mxy"; 
if curmodel=="xy"
    modeldir="xy_s";
    modelname="XY"; 
    sqrtN_vals=[16,32,64,128];
    L_vals=sqrtN_vals;
    T_str=".70"; 
    i_q_base=2;
elseif curmodel=="mxy"
    modeldir="mxy_3.00";
    modelname="MXY"; 
    sqrtN_vals=[16,32,64,128];
    L_vals=[18.5,37,74,148];
    T_str=".14";
    i_q_base=2;
end
T=str2double(T_str);
q_base = 2*pi*i_q_base/L_vals(1);

n_run=500;
n_avg=50;
n_samps=n_run/n_avg;

res_function = resolution_Gauss;

xmax=8;

fitstr="exp(-gamma*x/2)*(cos(omega_1*x)+gamma/2/omega_1*sin(omega_1*x))";
fitfunc=@(t,ga,om) exp(-ga*t/2).*(cos(om*t)+ga/2/om*sin(om*t));
fitfunc_DO_reciprocal=@(om,ga,om_0) 2*ga*om_0^2./((om.^2-om_0^2).^2+om.^2*ga^2);

n_trunc=10;
smoothweight=2;
smoothdepth=6;
% NF_spectrum=NF_Psi/q^(3-eta_cur);

%% 
for i_N=1:numel(sqrtN_vals)
    sqrtN=sqrtN_vals(i_N);
    fprintf("sqrtN = %d\n",sqrtN);
    filename="/data/scc/thobi/220201_ReducedSmapleStepDeltat/" + modeldir + ...
        "/sqrtN_" + sqrtN + "/T_" + T_str + "/samp_Dynamics.mat";
    S=load(filename, ...
        "averaging_times","qbin","absM_av","rbin", ...
        "gmperpmperp","gxx","gyy","gmparmpar");

    t=S.averaging_times;
    qbin=S.qbin;
    n_t=numel(t);
    n_q=numel(qbin);
    eta_cur=-2*log(S.absM_av)/log(sqrt(2)*sqrtN);
    absM_vals(i_N)=S.absM_av;
    
    [minval,argmin]=min(abs(qbin-q_base));
    i_q=argmin;
    q=qbin(i_q);
    
    gmperp=real(S.gmperpmperp)/sqrtN^2;
    gm=real(S.gxx + S.gyy)/sqrtN^2;
    gmpar=real(S.gmparmpar)/sqrtN^2;

    q_vals(i_N)=q;
    gmperp_cell{i_N}=gmperp;
    gm_cell{i_N}=gm;
    gmpar_cell{i_N}=gmpar;
    eta_vals(i_N)=eta_cur;

% eta_cur=-2*log(S.absM_av)/log(sqrt(2)*128);
% y_vals=sort([logspace(-2,1,500),.9955:.001:1.01]);
% NF_Psi=NelsonFisher_Psi(y_vals,eta_cur,n_trunc);



    cf=real(gmperp(i_q:n_q:end));
    fitob=fit(t(:),cf(:)/cf(1),fitstr,'StartPoint',[.1,.1]);
    cf_perp_cell{i_N}=cf;
    om_1_perp_vals(i_N)=abs(fitob.omega_1);
    ga_perp_vals(i_N)=fitob.gamma;
    c_vals(i_N) = fitob.omega_1 /q;
    
    cf_m=real(gm(i_q:n_q:end));
    fitob=fit(t(:),cf_m(:)/cf_m(1),fitstr,'StartPoint',[ga_perp_vals(i_N),om_1_perp_vals(i_N)]);
    cf_m_cell{i_N}=cf_m;
    om_1_m_vals(i_N)=abs(fitob.omega_1);
    ga_m_vals(i_N)=fitob.gamma;
    c_m_vals(i_N) = fitob.omega_1 /q;
    
    cf_par=real(gmpar(i_q:n_q:end));
    fitob=fit(t(:),cf_par(:)/cf_par(1),fitstr,'StartPoint',[.1,.1]);
    cf_par_cell{i_N}=cf_par;
    om_1_par_vals(i_N)=abs(fitob.omega_1);
    ga_par_vals(i_N)=fitob.gamma;
    c_par_vals(i_N) = fitob.omega_1 /q;

    tau=4/ga_perp_vals(i_N);
    res_vals=res_function(t,tau);
    [ft_vals,om_vals]=FT_correlation(t, cf .* res_vals /cf(1), 1e6);
    ft_perp_vals_cell{i_N}=real(ft_vals);
    om_perp_vals_cell{i_N}=om_vals;
    [ft_vals,om_vals]=FT_correlation(t, cf .* res_vals /cf(1), 1e6);

    [ft_m_vals,om_vals]=FT_correlation(t, cf_m .* res_vals /cf_m(1), 1e6);
    ft_m_vals_cell{i_N}=real(ft_m_vals);

    tau=4/ga_par_vals(i_N);
    res_vals=res_function(t,tau);
    [ft_par_vals,om_vals]=FT_correlation(t, cf_par .* res_vals /cf_par(1), 1e6);
    ft_par_vals_cell{i_N}=real(ft_par_vals);


end
y_vals=sort([logspace(-2,1,500),.9955:.001:1.01]);
NF_Psi=NelsonFisher_Psi(y_vals,eta_cur,n_trunc);
ind_om_pos=find(om_vals>=0);

%% perp FT
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
% c_map=linspecer(numel(sqrtN_vals));
c_map=colormap_sqrtN();
for i_N=1:numel(sqrtN_vals)
    sqrtN = sqrtN_vals(i_N);
    cf_cur=cf_perp_cell{i_N};
    ft_cur=ft_perp_vals_cell{i_N}(ind_om_pos) * cf_cur(1);
    om_vals_cur=om_vals(ind_om_pos);
    om_cur=om_1_perp_vals(i_N);
    dispname=sprintf("$N = (%d)^2$",sqrtN);
%     ft_real=real(ft_vals*cf(1)/sqrtN^2);
    ft_smooth=smoothen(abs(ft_cur),smoothweight,smoothdepth);
%     loglog(om_vals/om_cur,abs(ft_cur), ...
    loglog(om_vals_cur/om_cur,ft_smooth, ...
        '-','Displayname',dispname, ...
        'Color',c_map(i_N,:));
    hold on;
end
DO_vals=cf_cur(1)*fitfunc_DO_reciprocal(om_vals_cur,ga_perp_vals(end),om_1_perp_vals(end));
loglog(om_vals_cur/om_cur,DO_vals,'-.k','DisplayName','hydrodynamic fit', ...
    'LineWidth',1.2,...
    'HandleVisibility','off');

NF_cur=NF_Psi*ft_smooth(1)/NF_Psi(1);
loglog(y_vals,NF_cur,'--k','DisplayName','NF law (scaled)', ...
    'LineWidth',1.2,...
    'HandleVisibility','off');
% loglog(y_vals,NF_Psi/q^(3-eta_cur),'--k','DisplayName','NF law');
[minval, i_q_right] = min(abs(y_vals-xmax));
% ylim([.1*NF_Psi(i_q_right)/q^(3-eta_cur) max(NF_Psi/q^(3-eta_cur))]);
xlim_vec=[1e-1 xmax];
ylim_vec=[.1*NF_cur(i_q_right) 10*max(ft_cur)];
xlim(xlim_vec);
ylim(ylim_vec);

str="HD";
xval=4;
yval=cf_cur(1)*fitfunc_DO_reciprocal(xval*om_cur,ga_perp_vals(end),om_1_perp_vals(end));
xvec = [.6*xval, .9*xval];
yvec = [yval, yval];
xvec = log10(xvec/xlim_vec(1)) ./ log10(xlim_vec(2)/xlim_vec(1));
yvec = log10(yvec/ylim_vec(1)) ./ log10(ylim_vec(2)/ylim_vec(1));
valign='middle';
halign='right';
arrow_annot_HD=annotation('textarrow',xvec,yvec,'String',str,...
    'Units','normalized',...
    'LineWidth',.5,...
    'VerticalAlignment',valign, 'HorizontalAlignment',halign,...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_annotation,...
    'HeadStyle','plain','HeadWidth',2,'HeadLength',3);

str="NF";
xval=4;
[~,i_min]= min(abs(y_vals-xval));
xval=y_vals(i_min);
yval=NF_cur(i_min);
xvec = [xval, xval];
yvec = [6*yval, 1.1*yval];
xvec = log10(xvec/xlim_vec(1)) ./ log10(xlim_vec(2)/xlim_vec(1));
yvec = log10(yvec/ylim_vec(1)) ./ log10(ylim_vec(2)/ylim_vec(1));
valign='bottom';
halign='center';
arrow_annot_HD=annotation('textarrow',xvec,yvec,'String',str,...
    'Units','normalized',...
    'LineWidth',.5,...
    'VerticalAlignment',valign, 'HorizontalAlignment',halign,...
    'interpreter','latex',...
    'Color','black','FontSize', fontsize_annotation,...
    'HeadStyle','plain','HeadWidth',2,'HeadLength',3);

set(gca,'FontName', 'cmr12','FontSize', fontsize_annotation);
hlegend=legend('Location','SouthWest','Interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_annotation);
hXLabel = xlabel('$\omega / cq$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$S_{m\perp}(q,\omega)$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
h_axis = gca;


annotation_str = {modelname,sprintf('$T = %.3f$',T),sprintf('$q = %.3f$',q)};
dim=[.46 .16 .1 .2]; % Right of bottom legend
dim=[.16 .73 .1 .2]; % Above left bottom legend
h_annotation = annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');
figname=sprintf('%s/plots/%s_Smperp_T_%.3f_q_%.3f',basedir,curmodel,T,q);
fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

%% M FT
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
c_map=colormap_sqrtN();
for i_N=1:numel(sqrtN_vals)
    sqrtN = sqrtN_vals(i_N);
    cf_cur=cf_m_cell{i_N};
    ft_cur=ft_m_vals_cell{i_N}(ind_om_pos) * cf_cur(1);
    om_vals_cur=om_vals(ind_om_pos);
    om_cur=om_1_perp_vals(i_N);
    dispname=sprintf("$N = (%d)^2$",sqrtN);
    ft_smooth=smoothen(abs(ft_cur),smoothweight,smoothdepth);
    loglog(om_vals_cur/om_cur,ft_smooth, ...
        '-','Displayname',dispname, ...
        'Color',c_map(i_N,:));
    hold on;
end
DO_vals=cf_cur(1)*fitfunc_DO_reciprocal(om_vals_cur,ga_m_vals(end),om_1_m_vals(end));
loglog(om_vals_cur/om_cur,DO_vals,'-.k','DisplayName','DO fit', ...
    'LineWidth',1.2);
NF_cur=NF_Psi*ft_cur(1)/NF_Psi(1);
loglog(y_vals,NF_cur,'--k','DisplayName','NF law (scaled)', ...
    'LineWidth',1.2);
% loglog(y_vals,NF_Psi/q^(3-eta_cur),'--k','DisplayName','NF law');
[minval, i_q_right] = min(abs(y_vals-xmax));
% ylim([.1*NF_Psi(i_q_right)/q^(3-eta_cur) max(NF_Psi/q^(3-eta_cur))]);
xlim([1e-1 xmax]);
ylim([.1*NF_cur(i_q_right) 10*max(ft_cur)]);


hlegend=legend('Location','SouthWest','Interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_annotation);

hXLabel = xlabel('$\omega / cq$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$S_{m}(q,\omega)$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
h_axis = gca;

annotation_str = {modelname,sprintf('$T = %.3f$',T),sprintf('$q = %.3f$',q)};
dim=[.46 .16 .1 .2]; % Right of bottom legend
dim=[.16 .73 .1 .2]; % Above left bottom legend
h_annotation = annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

figname=sprintf('%s/plots/%s_Sm_T_%.3f_q_%.3f',basedir,curmodel,T,q);
fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end










%% Parallel FT
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
c_map=colormap_sqrtN();
for i_N=1:numel(sqrtN_vals)
    sqrtN = sqrtN_vals(i_N);
    cf_cur=cf_par_cell{i_N};
    ft_cur=ft_par_vals_cell{i_N}(ind_om_pos) * cf_cur(1);
    om_vals_cur=om_vals(ind_om_pos);
    om_cur=om_1_perp_vals(i_N);
    dispname=sprintf("$N = (%d)^2$",sqrtN);
%     ft_real=real(ft_vals*cf(1)/sqrtN^2);
    ft_smooth=smoothen(abs(ft_cur),smoothweight,smoothdepth);
%     loglog(om_vals/om_cur,abs(ft_cur), ...
    loglog(om_vals_cur/om_cur,ft_smooth, ...
        '-','Displayname',dispname, ...
        'Color',c_map(i_N,:));
    hold on;
end
% NF_cur=NF_Psi*ft_cur(1)/NF_Psi(1);
% loglog(y_vals,NF_cur,'--k','DisplayName','NF law (scaled)', ...
%     'LineWidth',1.2);
% loglog(y_vals,NF_Psi/q^(3-eta_cur),'--k','DisplayName','NF law');
[minval, i_q_right] = min(abs(y_vals-xmax));
% ylim([.1*NF_Psi(i_q_right)/q^(3-eta_cur) max(NF_Psi/q^(3-eta_cur))]);
xlim([1e-1 xmax]);
ylim([.1*NF_cur(i_q_right) 10*max(ft_cur)]);

hlegend=legend('Location','SouthWest','Interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_annotation);
hXLabel = xlabel('$\omega / cq$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
hYLabel = ylabel('$S_{m\parallel}(q,\omega)$','interpreter','latex', ...
    'FontName', 'cmr12','FontSize', fontsize_ax_labels);
h_axis = gca;

annotation_str = {modelname,sprintf('$T = %.3f$',T),sprintf('$q = %.3f$',q)};
dim=[.46 .16 .1 .2]; % Right of bottom legend
dim=[.16 .73 .1 .2]; % Above left bottom legend
h_annotation = annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'LineWidth', .5, ...
    'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

figname=sprintf('%s/plots/%s_Smpar_T_%.3f_q_%.3f',basedir,curmodel,T,q);
fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end