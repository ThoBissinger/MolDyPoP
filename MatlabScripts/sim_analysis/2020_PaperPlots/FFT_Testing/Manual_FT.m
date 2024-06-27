%% 0 Initialization
clear;
close all;
run ../initialization_script

i_q = 1;
res_factor=30;
res_function = resolution_Gauss;

load('/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_128/T_.14/samp_Dynamics.mat',...
    'averaging_times','gmperpmperp','qbin');

%% 1 Extract standard info

q = qbin(i_q);
n_q = numel(qbin);
t = averaging_times;
dt = t(2)-t(1);
cf = real(gmperpmperp(i_q:n_q:end));
cf = cf / cf(1);
n_t = numel(t);
coeffs_omega_1_cur=fit_DampedOscillator_RealSpace(t,real(cf),10,1,'omega_1');
gamma_cur = coeffs_omega_1_cur(1);
omega_1_cur = coeffs_omega_1_cur(2);
omega_0_cur = sqrt(omega_1_cur^2 + .25 * gamma_cur^2);
c_cur = omega_0_cur / q;

tau=res_factor/omega_1_cur;
resolution_vals=res_function(t,tau);
                    
[ft_vals,om_vals]=FT_correlation(t, cf .* resolution_vals, 0);

%% 2 Plot Function
figure
plot(t,cf);
hold on;
plot(t,cf.*resolution_vals);



%% 3 Plot Standard FFT with some damping
figure
semilogy(om_vals,ft_vals);
xlim([0 3*omega_0_cur]);

%% 4 Manual FFT
om_vals_man = linspace(0,3*omega_0_cur,1e3);
ft_vals_man = zeros(size(om_vals_man));
n_om = numel(om_vals_man);

for i_om = 1:n_om
    om = om_vals_man(i_om);
    indices = 1:n_t;
    ft_vals_man(i_om) = 2*dt*sum(cos(om * t(indices)).* cf(indices) .* resolution_vals(indices));
end

%% 5 Plot manual FFT
figure
semilogy(om_vals_man,ft_vals_man);
xlim([0 3*omega_0_cur]);
