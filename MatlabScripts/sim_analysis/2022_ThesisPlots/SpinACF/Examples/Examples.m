clear
pathbase="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/SpinACF/Examples";
cd(pathbase)
run('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/initialization_script.m')
saveswitch=1;
basedir=pathbase;

sqrtN_vals = [16 32 64 128];
sampfilename="samp_Dynamics";
    
    
% model="mxy"; model="xy"; model="xy_s";
model="mxy";
if model=="mxy"
    curmodel="mxy";
    sqrtN = 128;
    basedir="/data/scc/thobi/211201_LongerTime/mxy_3.00";
    T_str = ".17";
    T = str2double(T_str);
    i_q_vals = [3,6,9];
    xlim_vec_qt=[0 600];
    xlim_vec_qo=[0 .1];
    ylim_vec_qo=[1e1 Inf];
    ylim_vec_qo_lin=[0 Inf];
    smoothdepth=2;
    
elseif model=="xy_s"
    curmodel="xy_s";
    sqrtN = 128;
    basedir="/data/scc/thobi/211201_LongerTime/xy_s";
    T_str = ".89";
    T = str2double(T_str);
    i_q_vals = [2,3,6];
    xlim_vec_qt=[0 1e3];
    xlim_vec_qo=[0 .23];
    ylim_vec_qo=[5e1 Inf];
    ylim_vec_qo_lin=[0 Inf];
    smoothdepth=2;


elseif model=="xy"
    modeltag="XY";
    model="xy";
    curmodel="xy";
    sqrtN_vals = [16 32 64 128];
    L_vals = sqrtN_vals;
%     T_str = [".10" ".20" ".30" ".40" ".50" ".60" ".70" ".80" ".90" "1.00" "1.10" "1.20" "1.30" "1.40" "1.50" "1.60" "1.70" "1.80" "1.90" "2.00" ];
    T_str = [".10" ".30" ".50" ".70" ".90" "1.10" "1.30" "1.40" "1.50" "1.60" "1.70" "1.90" "2.10"];
    T_dirs = sprintfc("T_%s",T_str);
    T_vals = str2double(T_str);
    
    dir_reduced_dt="/data/scc/thobi/220201_ReducedSmapleStepDeltat/xy";
    simdir_helicity = "220201_ReducedSmapleStepDeltat";
%     dir_reduced_dt="/data/scc/thobi/210715_LinearTimeSampling/xy";
%     simdir_helicity = "210715_LinearTimeSampling";
    modeldir="xy";

    xlim_vec=[0 2];
    ylim_sigma_vec=[0 Inf];
    ylim_omega=[0 Inf];
    xlim_omega_inset=[0 2];
    ylim_omega_inset=[0 Inf];
    YTick_Upsilon_vec=.9:.05:1.1;

end

%load(sprintf('%s/%s_dataset',basedir,model));

fitfunc_symm="exp(-gamma*x/2)*(cos(omega_1*x) + gamma/2/omega_1*sin(omega_1*x))";
fitfunc_asymm="exp(-gamma*x/2)*cos(omega_1*x)";
f_symm=@(t,om,ga) exp(-ga*t/2).*(cos(om*t)+ga/2/om*sin(om*t));

simfile=sprintf("%s/sqrtN_%d/T_%s/samp_Dynamics",basedir,sqrtN,T_str);

S=load(simfile,"averaging_times","gmperpmperp","gxx","gyy","qbin");

res_fct=@(t,tau) exp(-t.^2/2/tau^2);

t=S.averaging_times;
gmperpmperp=S.gmperpmperp;
qbin=S.qbin;
n_q = numel(qbin);
q_vals=S.qbin(i_q_vals);
for ind_q = 1:numel(i_q_vals)
    i_q = i_q_vals(ind_q);
    fprintf('q = %.3f, ',q_vals(ind_q));
    cf{ind_q}=real(gmperpmperp(i_q:numel(qbin):end))/sqrtN^2;
    chi_vals(ind_q)=cf{ind_q}(1);
    fitob_sym=fit(t(:),cf{ind_q}(:)/chi_vals(ind_q),fitfunc_symm,'StartPoint',[1e-2, 1e-2], 'Lower', [0, 0]);
    ga_vals(ind_q)=fitob_sym.gamma;
    om_1_vals(ind_q)=abs(fitob_sym.omega_1);

    tau_vals(ind_q) = log(500)./ga_vals(ind_q);
    res_vals{ind_q} = res_fct(t,tau_vals(ind_q));
    [ft_cur,om_cur]=FT_correlation(t, cf{ind_q} .* res_vals{ind_q} /cf{ind_q}(1), 0);
    ft{ind_q}=real(ft_cur);
    om_vals{ind_q}=om_cur;

end

fprintf('\n');

%% time domain
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
for ind_q = 1:numel(i_q_vals)
    i_q = i_q_vals(ind_q);
    dispname=sprintf('$q = %.3f$',q_vals(ind_q));
    plot(t,cf{ind_q},'-','DisplayName',dispname);
    hold on;
end
xlim(xlim_vec_qt);
xlabel("$t$","interpreter",'Latex','FontSize',10,'FontName','cmr14')
ylabel("$C_{m\perp m\perp}(q,t)$","interpreter",'Latex','FontSize',10,'FontName','cmr14')
hLegend = legend('Location', 'SouthEast','interpreter','latex',...
    'NumColumns',1,'FontSize',8,'FontName','cmr14');    
figname=sprintf('%s/%s_C_mm_qt_%d_%.2f',pathbase,model,sqrtN,T);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');


%% frequency domain
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
for ind_q = 1:numel(i_q_vals)
    i_q = i_q_vals(ind_q);
    dispname=sprintf('$q = %.3f$',q_vals(ind_q));
    y_vals=smoothen(chi_vals(ind_q)*abs(ft{ind_q}),2,smoothdepth);
    semilogy(om_vals{ind_q},y_vals,'-','DisplayName',dispname);
    hold on;
end
xlim(xlim_vec_qo);
ylim(ylim_vec_qo);
xlabel("$\omega$","interpreter",'Latex','FontSize',10,'FontName','cmr14')
ylabel("$S_{m\perp m\perp}(q,\omega)$","interpreter",'Latex','FontSize',10,'FontName','cmr14')
hLegend = legend('Location', 'NorthEast','interpreter','latex',...
    'NumColumns',1,'FontSize',8,'FontName','cmr14');    

figname=sprintf('%s/%s_S_mm_qo_%d_%.2f',pathbase,model,sqrtN,T);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');

set(gca,'Yscale','linear');
ylim(ylim_vec_qo_lin);
figname=sprintf('%s/%s_S_mm_qo_%d_%.2f_lin',pathbase,model,sqrtN,T);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');