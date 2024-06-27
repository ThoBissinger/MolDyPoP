%% 0 Initialization
clear
run ../initialization_script;
saveswitch=1;
basedir="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/NF_FFT";


modeldir="xy_s"; modelname="XY"; q=0.4048;
sqrtN_vals=128;
T_str=".85"; T=str2double(T_str);
filename="/data/scc/thobi/220201_ReducedSmapleStepDeltat/" + modeldir + ...
    "/sqrtN_" + sqrtN + "/T_" + T_str + "/samp_Dynamics.mat";

i_q=2^2*6;

n_run=500;
n_avg=50;
n_samps=n_run/n_avg;

S=load(filename, ...
    "averaging_times","qbin","absM_av","rbin", ...
    "gmperpmperp","gxx","gyy","gmparmpar");

t=S.averaging_times;
q_vals=S.qbin;
n_t=numel(t);
n_q=numel(q_vals);
eta_cur=-2*log(S.absM_av)/log(sqrt(2)*sqrtN);
q=q_vals(i_q);

gmperp=real(S.gmperpmperp);
gm=real(S.gxx + S.gyy);
gmpar=real(S.gmparmpar);


res_function = resolution_Gauss;

fitstr="exp(-gamma*x/2)*(cos(omega_1*x)+gamma/2/omega_1*sin(omega_1*x))";
fitfunc=@(t,ga,om) exp(-ga*t/2).*(cos(om*t)+ga/2/om*sin(om*t));

n_trunc=10;

y_vals=sort([logspace(-2,1,500),.9955:.001:1.01]);
NF_Psi=NelsonFisher_Psi(y_vals,eta_cur,n_trunc);
% NF_spectrum=NF_Psi/q^(3-eta_cur);

%% 

% eta_cur=-2*log(S.absM_av)/log(sqrt(2)*128);
% y_vals=sort([logspace(-2,1,500),.9955:.001:1.01]);
% NF_Psi=NelsonFisher_Psi(y_vals,eta_cur,n_trunc);



cf=real(gmperp(i_q:n_q:end));
fitob=fit(t(:),cf(:)/cf(1),fitstr,'StartPoint',[.1,.1]);
om=fitob.omega_1;
ga=fitob.gamma;
c = om /q;

cf_m=real(gm(i_q:n_q:end));
fitob=fit(t(:),cf_m(:)/cf_m(1),fitstr,'StartPoint',[.1,.1]);
om_m=fitob.omega_1;
ga_m=fitob.gamma;
c_m = om_m /q;

cf_par=real(gmpar(i_q:n_q:end));
fitob=fit(t(:),cf_par(:)/cf_par(1),fitstr,'StartPoint',[.1,.1]);
om_par=fitob.omega_1;
ga_par=fitob.gamma;
c_par = om_par /q;


figure
plot(t,[cf_m;cf;cf_par]);

%% perp FT
tau=1/ga;
res_vals=res_function(t,tau);
[ft_vals,om_vals]=FT_correlation(t, cf .* res_vals /cf(1), 1e6);

figure
ft_real=real(ft_vals*cf(1)/sqrtN^2);
ft_smooth=smoothen(abs(ft_real),1,3);
loglog(om_vals/om,ft_real,'-r','Displayname','FT');
hold on;
loglog(om_vals/om,ft_smooth,'-b','Displayname','FT smooth');
loglog(y_vals,NF_Psi/q^(3-eta_cur),'-k','DisplayName','NF law');
xlim([1e-1 6]);


%% M FT
tau=1/ga_m;
res_vals=res_function(t,tau);
[ft_m_vals,om_vals]=FT_correlation(t, cf_m .* res_vals /cf_m(1), 1e6);

figure
ft_m_real=real(ft_m_vals*cf_m(1)/sqrtN^2);
ft_m_smooth=smoothen(abs(ft_m_real),2,5);
loglog(om_vals/om,ft_m_real,'-r','Displayname','FT');
hold on;
loglog(om_vals/om,ft_m_smooth,'-b','Displayname','FT smooth');
loglog(y_vals,NF_Psi/q^(3-eta_cur),'--k','DisplayName','NF law');
xlim([1e-1 6]);
legend show;

%% Parallel FT
tau=1/ga_par;
res_vals=res_function(t,tau);
[ft_par_vals,om_vals]=FT_correlation(t, cf_par .* res_vals /cf_par(1), 1e6);

figure
ft_par_real=real(ft_par_vals*cf_par(1)/sqrtN^2);
ft_par_smooth=smoothen(abs(ft_par_real),2,5);
loglog(om_vals/om,ft_par_real,'-r','Displayname','FT');
hold on;
loglog(om_vals/om,ft_par_smooth,'-b','Displayname','FT smooth');
loglog(y_vals,NF_Psi/q^(3-eta_cur),'--k','DisplayName','NF law');
xlim([1e-1 6]);
legend show;