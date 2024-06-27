%% 0 Initialization
clear
run ../../initialization_script

filename_short='/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00/sqrtN_128/T_.11/samp_Dynamics.mat';
filename_long='/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_128/T_.11/samp_Dynamics.mat';

i_q=6;

load(filename_long,'qbin','averaging_times','gmperpmperp')
cf_long=real(gmperpmperp(i_q:numel(qbin):end));
cf_long = cf_long / cf_long(1);
t_long=averaging_times;

load(filename_short,'qbin','averaging_times','gmperpmperp')
cf=real(gmperpmperp(i_q:numel(qbin):end));
cf = cf / cf(1);
t=averaging_times;

[t_c,y_c] = continue_cf(t,cf,3);

%% 1 FFT stuff
tau=3000;
res_vals=resolution_Gauss(t,tau);
res_vals_c=resolution_Gauss(t_c,tau);
res_vals_long=resolution_Gauss(t_long,tau);


[ft_vals,om_vals]=FT_correlation(t, cf.*res_vals, 1e5);
[ft_vals_c,om_vals_c]=FT_correlation(t_c, y_c.*res_vals_c, 1e5);
[ft_vals_long,om_vals_long]=FT_correlation(t_long, cf_long.*res_vals_long, 1e6);

[ft_max,i_max] = max(ft_vals_long);
om_max = om_vals_long(i_max);

ft_smooth=smoothen(abs(real(ft_vals)),1,1);
ft_smooth_c=smoothen(abs(real(ft_vals_c)),1,1);
ft_smooth_long=smoothen(abs(real(ft_vals_long)),1,3);



%% 2 Simple Plot of cf and continued cf
figure
c_map = linspecer(3);
set(0,'DefaultLineLineWidth',2)
plot(t_long,cf_long,'DisplayName','Regular CF',...
    'Color',c_map(1,:));
hold on;
plot(t_c,y_c,'DisplayName','Continued CF',...
    'Color',c_map(2,:));
plot(t,cf,'DisplayName','Short-Time CF',...
    'Color',c_map(3,:));
hold off;
xlim([0 3000]);
legend show;

%% 3 Plot of FFT
figure
c_map = linspecer(3);
set(0,'DefaultLineLineWidth',2)
loglog(om_vals_long / om_max,ft_smooth_long,'DisplayName','Regular CF',...
    'Color',c_map(1,:));
hold on;
loglog(om_vals_c / om_max,ft_smooth_c,'DisplayName','Continued CF',...
    'Color',c_map(2,:));
loglog(om_vals / om_max,abs(real(ft_vals)),'DisplayName','Short-Time CF',...
    'Color',c_map(3,:));
hold off;
xlim([1e-1 20]);
legend show;