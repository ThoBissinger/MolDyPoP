% clear
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
if ispc
    fig_base = "C:/Users/tbiss/Thesis/MatlabScripts/sim_analysis/2020_PaperPlots";
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\ExternalCodes\')
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\sim_analysis\FitFuncs\')
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\sim_analysis\FourierTransform\')
%     addpath('C:\Users\tbiss\Thesis\MatlabScripts\sim_analysis\data_access\')
else
    fig_base = "/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots";
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes')
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/FitFuncs')
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/FourierTransform')
%     addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/data_access/')
end
% addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
% addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;


mxydata_name='mxy/rho_3.00_dynamics_LinearTime.mat';
xydata_name='xy/xy_dynamics_LinearTime.mat';
fmxydata_name='fmxy/fmxy_dynamics_LinearTime.mat';
xysdata_name='xy/eq_xy_s.mat';

mxydata_AdjustedTime_name='mxy/mxy_dynamics_AdjustedTime.mat';
xydata_AdjustedTime_name='xy/xy_dynamics_AdjustedTime.mat';
fmxydata_AdjustedTime_name='fmxy/fmxy_dynamics_AdjustedTime.mat';
xysdata_AdjustedTime_name='xy/xy_s_dynamics_AdjustedTime.mat';

mxyfit_TCF_q_name='mxy/rho_3.00_TCF_q_fit_LinearTime.mat';
xyfit_TCF_q_name='xy/xy_TCF_q_fit_LinearTime.mat';
fmxyfit_TCF_q_name='fmxy/fmxy_TCF_q_fit_LinearTime.mat';
xysfit_TCF_q_name='xy/xy_s_TCF_q_fit_mperp_AdjustedTime';

mxyfit_name='mxy/rho_3.00_CritExpFit.mat';
xyfit_name='xy/xy_CritExpFit.mat';
fmxyfit_name='fmxy/fmxy_CritExpFit.mat';

mxyfit_FS_name='mxy/rho_3.00_DataCollapse_SCF.mat';
xyfit_FS_name='xy/xy_DataCollapse_SCF.mat';
fmxyfit_FS_name='fmxy/fmxy_DataCollapse_SCF.mat';

mxyfit_FSMag_name='mxy/rho_3.00_etafit_FSMag.mat';
xyfit_FSMag_name='xy/xy_etafit_FSMag.mat';
fmxyfit_FSMag_name='fmxy/fmxy_etafit_FSMag.mat';
xysfit_FSMag_name='xy/xy_s_etafit_FSMag.mat';

mxydata_LepriRuffo_name='mxy/rho_3.00_LepriRuffo_extended.mat';
xydata_LepriRuffo_name='xy/xy_LepriRuffo_extended.mat';
fmxydata_LepriRuffo_name='fmxy/fmxy_LepriRuffo_extended.mat';

mxydata_LepriRuffo_extended_name='mxy/rho_3.00_LepriRuffo_extended.mat';
xydata_LepriRuffo_extended_name='xy/xy_LepriRuffo_extended.mat'; % doesn't exist yet
fmxydata_LepriRuffo_extended_name='fmxy/fmxy_LepriRuffo_extended.mat'; % doesn't exist yet

resolution_Laplace=@(t,tau) exp(-t/tau);
resolution_Gauss=@(t,sigma) exp(-.5*t.^2/sigma^2);
resolution_Quartic=@(t,sigma) exp(-(t/sigma).^4/(24));
resolution_Laplace_pleateau=@(t,tau,t_max) (t <= t_max) + (t > t_max) .*exp(-(t-t_max)/tau);

fitfunc_DO=@(times,corrfunc_1,c) corrfunc_1 * exp(-c(1) * times/2) .* (cos(c(2) * times) + .5*c(1)/c(2)*sin(c(2) * times));  
fitfunc_DO_reciprocal=@(om,corrfunc_1,c) 2*corrfunc_1*c(1)*c(2)^2 ./ ...
    ( ( om.^2 - c(2)^2 ).^2 + om .^2 * c(1)^2 );
om_1=@(omega_0,gamma) sqrt(omega_0^2 - .25*gamma^2);
om_2=@(omega_0,gamma) sqrt(omega_0^2 - .5*gamma^2);
    
fontsize_annotation = 7; % For annotations within the figure, including legend
fontsize_labels = 10; % For larger labels like titles etc
fontsize_axis = 7; % For the axes, that is the ticks.
fontsize_ax_labels = 10; % For axis labels, xlabel and ylabel
fontsize_titles = 12; % For image titles
fontsize_subfiglabels = 7; % For the labels of subfigures, like (a), (b), etc

columnwidth_cm = 9;
oneandahalfpagewidth_cm = 14;
pagewidth_cm = 19;

pageheight_cm = 24.7;

dpi_png = 300; % Should be ramped up for publication-worthy stuff

n_period = 10; weightexp = 1; % For fitting fit_DampedOscillator_RealSpace(t,cf,n_period,weightexp,'omega_1');
