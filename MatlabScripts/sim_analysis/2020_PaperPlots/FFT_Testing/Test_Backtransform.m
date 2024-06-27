close all;
% clear all;


if ispc
    fig_base = "C:/Users/tbiss/Thesis/MatlabScripts/sim_analysis/2020_PaperPlots/FFT_Testing/plots/Test_Backtransform";
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\ExternalCodes\')
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\sim_analysis\FitFuncs\')
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\sim_analysis\FourierTransform\')
else
    fig_base = "/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/FFT_Testing/plots/Test_Backtransform";
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScriptsExternalCodes')
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScriptssim_analysis/FitFuncs')
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScriptssim_analysis/FourierTransform')
end
saveswitch = 1;


DO_func=@(t,omega,gamma) exp(-.5*gamma * t) .* (cos(omega*t) + .5*gamma/omega*sin(omega*t));
% FFT_exact=@(om,omega_0,gamma) 3*gamma*omega_0^2 ./ ...
%     ( ( om.^2 - omega_0^2 ).^2 + om .^2 * 12*gamma^2 );
FFT_exact=@(om,omega_0,gamma) 2*gamma*omega_0^2 ./ ...
    ( ( om.^2 - omega_0^2 ).^2 + om .^2 * gamma^2 );
resolution_func=@(t,sigma) exp(-.5 * t.^2 /sigma^2);
om_1=@(omega_0,gamma) sqrt(omega_0^2 - .25*gamma^2);
om_2=@(omega_0,gamma) sqrt(omega_0^2 - .5*gamma^2);
cos_func=@(t,omega_0,sigma) exp(-.5 * t.^2/sigma^2) .* (cos(omega_0*t));

%% 1 Damped oscillator in Fourier space
figure(1);
dt=.1;
om_vals=linspace(-pi/dt,pi/dt,1e6);
omega=1e-2;
gamma_vals = [.5, sqrt(2),2,2.5]*omega;
labels=["$\omega_1^2 > 0,\ \omega_2^2>0$", "$\omega_1^2 > 0,\ \omega_2^2 = 0$", "$\omega_1^2 = 0,\ \omega_2^2 < 0$", "$\omega_1^2 < 0,\ \omega_2^2 < 0$"];
c_map=linspecer(numel(gamma_vals));
for i = 1 :numel(gamma_vals)
    gamma = gamma_vals(i);
    fprintf('om_1^2 = %.6f,    om_2^2 = %.6f\n', om_1(omega,gamma)^2, om_2(omega,gamma)^2);
    fvals=FFT_exact(om_vals,omega,gamma);
    plot(om_vals,fvals,...
        'DisplayName',labels(i),...
        'LineWidth',2,...
        'Color',c_map(i,:));
    hold on;
end
xlim([-.04,.04]);
legend show
set(legend,'interpreter','latex')
xlabel('$\omega$','interpreter','latex');
ylabel('$\hat{h}(\omega)$','interpreter','latex');

figname=sprintf('%s/DO_omega_space_examples',fig_base);
fprintf('Completed figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end


%% 2 Back transform into real space with fftshift
figure(2)
gamma_vals = [.5, sqrt(2),2,2.5]*omega;
labels=["$\omega_1^2 > 0,\ \omega_2^2>0$", "$\omega_1^2 > 0,\ \omega_2^2 = 0$", "$\omega_1^2 = 0,\ \omega_2^2 < 0$", "$\omega_1^2 < 0,\ \omega_2^2 < 0$"];
c_map=linspecer(numel(gamma_vals));
for i = 1 :numel(gamma_vals)
    gamma = gamma_vals(i);
    fvals=FFT_exact(om_vals,omega,gamma);
    [realfunc,t_vals]=FT_correlation(om_vals+min(om_vals),fftshift(fvals),0);
    plot(t_vals,real(realfunc)/4/pi,...
        'DisplayName',labels(i),...
        'LineWidth',2,...
        'Color',c_map(i,:));
%         'Displayname','$\hat{\hat{h}}(t)$');
    % plot(t_vals,real(realfunc)/4/pi.*((-1).^(0:numel(t_vals)-1)));
    hold on;
    omega_1 = om_1(omega,gamma);
    if (omega_1 == 0)
       omega_1 = 1e-12;
    end
    plot(t_vals(1:1e3:end),DO_func(abs(t_vals(1:1e3:end)),omega_1,gamma),...
        'DisplayName',labels(i),'HandleVisibility','off',...
        'LineWidth',2,...
        'LineStyle','none',...
        'Marker','^',...
        'Color',c_map(i,:));
%         'Displayname','$h(t)$');
end
xlim([0,2e3]);
hold off;
legend show
set(legend,'interpreter','latex');
xlabel('$t$','interpreter','latex');
ylabel('$h(t)$','interpreter','latex');

figname=sprintf('%s/DO_t_space_examples',fig_base);
fprintf('Completed figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end


%% 3 Option 1 of double transform, removing the negative values
% This quarters the number of data points, halving the window as well as
% doubling the minimal resolution. 
figure(3)

fvals=cos_func(t_vals,omega,sigma);
t_indices=find(t_vals>=0);
[ft_vals,om_vals]=FT_correlation(t_vals(t_indices),fvals(t_indices), 0);
om_indices=find(om_vals>=0);
[fvals_fft,t_vals_fft]=FT_correlation(om_vals(om_indices),real(ft_vals(om_indices)), 0);

plot(t_vals,fvals);
hold on;
plot(t_vals_fft,real(fvals_fft)/2/pi);
hold off;
fprintf('Double FFT with half-sided evaluation.  numel(t_vals) = %d,  numel(t_vals_fft) = %d\n', numel(t_vals), numel(t_vals_fft)); 

%% 4 Option 2 of double transform using fftshift
% Careful, this introduces a factor of 2 for each transform since the
% periodic continuation of the function leads to a double coverage of the
% interval in question. However, it leads to no loss in the number of
% time-domain points, while the other reduces the number of points
figure(4)
fvals=cos_func(t_vals,omega,sigma);
[ft_vals,om_vals]=FT_correlation(t_vals+min(t_vals),fftshift(fvals), 0);
[fvals_fft,t_vals_fft]=FT_correlation(om_vals+min(om_vals),fftshift(real(ft_vals)), 0);


plot(t_vals,fvals);
hold on;
plot(t_vals_fft,real(fvals_fft)/2/pi/4);
hold off;

fprintf('Double FFT with fftshift function.  numel(t_vals) = %d,  numel(t_vals_fft) = %d\n', numel(t_vals), numel(t_vals_fft)); 

figure(2)