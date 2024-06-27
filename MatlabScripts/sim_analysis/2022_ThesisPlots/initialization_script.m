% clear
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
if ispc
    fig_base = "C:/Users/tbiss/Thesis/MatlabScripts/sim_analysis/2022_ThesisPlots";
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\ExternalCodes\')
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\sim_analysis\FitFuncs\')
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\sim_analysis\FourierTransform\')
%     addpath('C:\Users\tbiss\Thesis\MatlabScripts\sim_analysis\data_access\')
else
    fig_base = "/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots";
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes')
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/FitFuncs')
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/FourierTransform')
%     addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/data_access/')
end
% addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));

% resolution_Laplace=@(t,tau) exp(-t/tau);
% resolution_Gauss=@(t,sigma) exp(-.5*t.^2/sigma^2);
% resolution_Quartic=@(t,sigma) exp(-(t/sigma).^4/(24));
% resolution_Laplace_pleateau=@(t,tau,t_max) (t <= t_max) + (t > t_max) .*exp(-(t-t_max)/tau);
% 
% fitfunc_DO=@(times,corrfunc_1,c) corrfunc_1 * exp(-c(1) * times/2) .* (cos(c(2) * times) + .5*c(1)/c(2)*sin(c(2) * times));  
% fitfunc_DO_reciprocal=@(om,corrfunc_1,c) 2*corrfunc_1*c(1)*c(2)^2 ./ ...
%     ( ( om.^2 - c(2)^2 ).^2 + om .^2 * c(1)^2 );
% om_1=@(omega_0,gamma) sqrt(omega_0^2 - .25*gamma^2);
% om_2=@(omega_0,gamma) sqrt(omega_0^2 - .5*gamma^2);
    
fontsize_annotation = 7; % For annotations within the figure, including legend
fontsize_labels = 10; % For larger labels like titles etc
fontsize_axis = 7; % For the axes, that is the ticks.
fontsize_ax_labels = 10; % For axis labels, xlabel and ylabel
fontsize_titles = 12; % For image titles
fontsize_subfiglabels = 7; % For the labels of subfigures, like (a), (b), etc

% fontsize_annotation = 9; % For annotations within the figure, including legend
% fontsize_labels = 12; % For larger labels like titles etc
% fontsize_axis = 9; % For the axes, that is the ticks.
% fontsize_ax_labels = 12; % For axis labels, xlabel and ylabel
% fontsize_titles = 12; % For image titles
% fontsize_subfiglabels = 9; % For the labels of subfigures, like (a), (b), etc

columnwidth_cm = 9;
oneandahalfpagewidth_cm = 14;
pagewidth_cm = 19;

pageheight_cm = 24.7;

dpi_png = 300; % Should be ramped up for publication-worthy stuff