%% 0 Initialization
run ../initialization_script
basedir=sprintf('%s/NF_Testing/Test_Matrix_FFT',fig_base);

%% 1 Simple Example
close all;
n_points=1e3;
tmax=1e2;
t=linspace(0,tmax,n_points);
dt=t(2)-t(1);
om_vals = linspace(-pi/dt,pi/dt,n_points);
res_vals=resolution_Gauss(t,tmax/10);
u=[cos(.2*t).*res_vals;cos(.5*t).*res_vals];
u_om = 2*dt*real(fftshift(fft(u'),1))';
i_u = 2;
u_om_exact = 2*dt*real(fftshift(fft(u(i_u,:))));
% u_om_exact = 2*dt*real(fft(u(i_u,:)));
plot(om_vals,u_om_exact - u_om(i_u,:));
% hold on;
% plot(om_vals,u_om(i_u,:));

% xlim([0 1]);
