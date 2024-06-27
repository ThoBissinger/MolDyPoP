close all;
% clear all;


if ispc
    fig_base = "C:/Users/tbiss/Thesis/MatlabScripts/sim_analysis/2020_PaperPlots/FFT_Testing/plots/FFT_DO_FFT_test";
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\ExternalCodes\')
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\sim_analysis\FitFuncs\')
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\sim_analysis\FourierTransform\')
else
    fig_base = "/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/FFT_Testing/plots/FFT_DO_FFT_test";
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScriptsExternalCodes')
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScriptssim_analysis/FitFuncs')
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScriptssim_analysis/FourierTransform')
end
saveswitch = 0;

tmax=1e4;


DO_func=@(t,omega,gamma) exp(-.5*gamma * t) .* (cos(omega*t) + .5*gamma/omega*sin(omega*t));
% FFT_exact=@(om,omega_0,gamma) 3*gamma*omega_0^2 ./ ...
%     ( ( om.^2 - omega_0^2 ).^2 + om .^2 * 12*gamma^2 );
FFT_exact=@(om,omega_0,gamma) 2*gamma*omega_0^2 ./ ...
    ( ( om.^2 - omega_0^2 ).^2 + om .^2 * gamma^2 );
resolution_func=@(t,sigma) exp(-.5 * t.^2 /sigma^2);
om_1=@(omega_0,gamma) sqrt(omega_0^2 - .25*gamma^2);

%% 1 Examples for Damped Oscillators
figure(1);

dt=1e1;
dt_str='$\Delta t = 10$';
times=0:dt:tmax;
tmax=1e4;
tmax_str='$t_{\textrm{max}} = 10^4$';

omega_vals=[2e-2,1e-2,4e-2];
gamma_vals=[5e-4,2e-4,3e-4];
names={'$\omega_0 = 2 \cdot 10^{-2}$, $\gamma = 5 \cdot 10^{-4}$',...
    '$\omega_0 = 1 \cdot 10^{-2}$, $\gamma = 2 \cdot 10^{-4}$',...
    '$\omega_0 = 4 \cdot 10^{-2}$, $\gamma = 3 \cdot 10^{-4}$'};
c_map=linspecer(numel(omega_vals));
for i = 1:numel(omega_vals)
    omega=omega_vals(i);
    gamma=gamma_vals(i);
    dispname=sprintf('$\\omega_0 = %d \\cdot 10^{-2},\\ \\gamma = %d \\cdot 10^{-4}$',omega*1e2,gamma*1e4);
    plot(times,DO_func(times,om_1(omega,gamma),gamma)-2.5*(i-1),...
        'DisplayName',names{i}, ...
        'Color',c_map(i,:), ...
        'LineWidth',2);
    hold on;
end
legend show
set(legend,'interpreter','latex','location','north')
set(gca,'YTick',[])
ylim([-7, 3]);
xlabel('t');
ylabel('$h(t)$','Interpreter','latex');

pos = get(gca, 'position');
dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)];
str = {dt_str,tmax_str};
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'interpreter','latex');

figname=sprintf('%s/DO_Examples',fig_base);
fprintf('Completed figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

% return


%% 2 FFT without noise
figure(2);

dt=1e1;
dt_str='$\Delta t = 10$';
times=0:dt:tmax;
tmax=1e4;
tmax_str='$t_{\textrm{max}} = 10^4$';


omega_vals=[2e-2,1e-2,4e-2];
gamma_vals=[5e-4,2e-4,3e-4];
names={'$\omega_0 = 2 \cdot 10^{-2}$, $\gamma = 5 \cdot 10^{-4}$',...
    '$\omega_0 = 1 \cdot 10^{-2}$, $\gamma = 2 \cdot 10^{-4}$',...
    '$\omega_0 = 4 \cdot 10^{-2}$, $\gamma = 3 \cdot 10^{-4}$'};

c_map=linspecer(numel(omega_vals));
for i = 1:numel(omega_vals)
    omega=omega_vals(i);
    gamma=gamma_vals(i);
    [ft_vals,om_vals]=FT_correlation(times,DO_func(times,om_1(omega,gamma),gamma), 0);
    dispname=sprintf('$\\omega_0 = %d \\cdot 10^{-2},\\ \\gamma = %d \\cdot 10^{-4}$',omega*1e2,gamma*1e4);
    plot(om_vals,real(ft_vals),...
        'DisplayName',names{i}, ...
        'Color',c_map(i,:), ...
        'LineWidth',2);
    hold on;
end
legend show
set(legend,'interpreter','latex','location','north')
xlabel('$\omega$','Interpreter','latex');
set(gca,'YScale','lin')
xlim([-.07,.07]);
ylim([-5e2,8e3]);

pos = get(gca, 'position');
dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)];
str = {dt_str,tmax_str};
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'interpreter','latex');
% 
figname=sprintf('%s/DO_FFT',fig_base);
fprintf('Completed figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end




%% 3 FFT Sigma compare
figure(3);
dt=1e1;
dt_str='$\Delta t = 10$';
times=0:dt:tmax;

omega=6e-2;
omega_str='$\omega_0 = 6 \cdot 10^{-2}$';
tmax=1e4;
tmax_str='$t_{\textrm{max}} = 10^4$';

omega_vals=[2e-2,1e-2,4e-2];
gamma_vals=[5e-4,2e-4,1e-4];
omega_str='$\omega_0 = 4 \cdot 10^{-2}$';
gamma_str='$\gamma = 10^{-4}$';

sigma_vals=[1e3,3e3,5e3,10e3,20e4];
sig_str={'$\sigma = 1\cdot 10^{3}$', '$\sigma = 3\cdot 10^{3}$', '$\sigma = 5\cdot 10^{3}$', '$\sigma = 1\cdot 10^{4}$', '$\sigma = 1\cdot 10^{5}$'};
rand_strength=0;
c_map=linspecer(numel(sigma_vals));
for i = 1:numel(sigma_vals)
    sigma=sigma_vals(i);
    for j = 1 :numel(omega_vals)
        if j == 1
            vis_switch = 'on';
        else
            vis_switch = 'off';
        end
        omega=omega_vals(j);
        gamma=gamma_vals(j);
        [ft_vals,om_vals]=FT_correlation(times,DO_func(times,om_1(omega,gamma),gamma).*resolution_func(times,sigma), 0);
    %     dispname=sprintf('$\\omega = %.1f \\cdot 10^{-2}, \\sigma = %.1f \\cdot 10^{3}$',omega*1e2,sigma/1e3);
%         om_select=find(( abs(real(ft_vals))>20) .* (om_vals > 0));
%         plot(om_vals(om_select),real(ft_vals(om_select)),...
        plot(om_vals,real(ft_vals),...
            'DisplayName',sig_str{i}, ...
            'Color',c_map(i,:), ...
            'HandleVisibility',vis_switch, ...
            'LineWidth',2,'LineStyle','-');
        hold on;
    end
end
xlim([.03,.05]);
legend show
set(legend,'interpreter','latex','location','east')
xlabel('$\omega$','Interpreter','latex');
ylabel('$\hat{h}(\omega)$','Interpreter','latex');
set(gca,'YScale','lin')

pos = get(gca, 'position');
dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)];
str = {omega_str,gamma_str,dt_str,tmax_str};
% str = {'$\omega_0 = 6 \cdot 10^{-2}$',sprintf('$\\Delta t = %.1f$',dt),'$t_{\textrm{max}} = 10^4$'};
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'interpreter','latex');
% 
figname=sprintf('%s/DO_FFT_sigma_compare',fig_base);
fprintf('Completed figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end


%% 4 FFT with noisy data
figure(4);
tmax=1e4;
dt=1e0;
dt_str='$\Delta t = 1$';
omega=6e-2;
omega_str='$\omega_0 = 6 \cdot 10^{-2}$';
gamma=5e-4;
gamma_str='$\gamma = 5\cdot 10^{-4}$';
sigma=3e3;
sigma_str='$\sigma = 3\cdot 10^{3}$';
tmax=1e4;
tmax_str='$t_{\textrm{max}} = 10^4$';
times=0:dt:tmax;

rand_strength=[0,0.1,.8];
c_map=linspecer(numel(rand_strength));
for i = 1:numel(rand_strength)
    omega=omega_vals(1);
    sigma=sigma_vals(1);
    [ft_vals,om_vals]=FT_correlation(times,DO_func(times,om_1(omega,gamma),gamma).*resolution_func(times,sigma)...
        + rand_strength(i)*rand(size(times)), 0);
    dispname=sprintf('$\\xi = %.1f$',rand_strength(i));
    if rand_strength(i) == 0
        lw = 2.5;
    else
        lw = 1.2;
    end
    plot(om_vals,abs(real(ft_vals)),...
        'DisplayName',dispname, ...
        'Color',c_map(i,:), ...
        'LineWidth',lw);
    hold on;
end
xlim([0,.15]);
ylim([1e-1,1e4]);
legend show
set(legend,'interpreter','latex','location','northwest')
xlabel('$\omega$','Interpreter','latex');
ylabel('$|\hat{h}(\omega)|$','Interpreter','latex');
set(gca,'YScale','log')

pos = get(gca, 'position');
dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)];
str = {omega_str,gamma_str,sigma_str,dt_str,tmax_str};
% str = {'$\omega_0 = 6 \cdot 10^{-2}$','$\sigma = 5\cdot 10^{-3}$',sprintf('$\\Delta t = %.1f$',dt),'$t_{\textrm{max}} = 10^4$'};
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'interpreter','latex');

% 
figname=sprintf('%s/DO_FFT_noisy',fig_base);
fprintf('Completed figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end


%% 5 Effect of dt
figure(5);
omega=4e-2;
omega_str='$\omega_0 = 4 \cdot 10^{-2}$';
gamma=1e-3;
gamma_str='$\gamma = 1\cdot 10^{-3}$';
sigma=3e3;
sigma_str='$\sigma = 3\cdot 10^{3}$';
tmax=1e4;
tmax_str='$t_{\textrm{max}} = 10^4$';

dt_vals=[1e1,1e0,1e-1,1e-2,1e-3];
c_map=linspecer(numel(dt_vals)+1);
for i = 1:numel(dt_vals)
    dt=dt_vals(i);
    times=0:dt:tmax;
    [ft_vals,om_vals]=FT_correlation(times,DO_func(times,om_1(omega,gamma),gamma).*resolution_func(times,sigma), 0);
    dispname=sprintf('$\\Delta t = %.3f$',dt);
    plot(om_vals,abs(FFT_exact(om_vals,omega,gamma)-real(ft_vals)),...
        'DisplayName',dispname, ...
        'Color',c_map(i,:), ...
        'LineWidth',2);
    hold on;
end
% plot(om_vals,FFT_exact(om_vals,omega,gamma),...
%         'DisplayName','exact', ...
%         'Color',c_map(end,:), ...
%         'LineWidth',2);
xlim([-.3,.3]);
ylim([1e-4,1e4]);
xlabel('$\omega$','Interpreter','latex');
ylabel('$|\hat{h}(\omega)|$','Interpreter','latex');
set(gca,'YScale','log')
legend show
set(legend,'interpreter','latex','location','northwest')

pos = get(gca, 'position');
dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)];
str = {omega_str,gamma_str,sigma_str,tmax_str};
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'interpreter','latex');
% 
figname=sprintf('%s/DO_FFT_dt',fig_base);
fprintf('Completed figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end



%% 6 Effect of spline
figure(6);
omega=4e-2;
omega_str='$\omega_0 = 4 \cdot 10^{-2}$';
gamma=1e-3;
gamma_str='$\gamma = 1\cdot 10^{-3}$';
sigma=3e3;
sigma_str='$\sigma = 3\cdot 10^{3}$';
tmax=1e4;
tmax_str='$t_{\textrm{max}} = 10^4$';
dt=10;
dt_str='$\Delta t = 10$';
times=0:dt:tmax;

n_points_vals=[0,1e4,1e5,5e5,1e6];
str_points={'no spline', '$n_{\textrm{spline}} = 10^4$', '$n_{\textrm{spline}} = 10^5$', '$n_{\textrm{spline}} = 5 \cdot 10^5$','$n_{\textrm{spline}} = 10^6$'};
c_map=linspecer(numel(n_points_vals)+1);
for i = 1:numel(n_points_vals)
    n_points=n_points_vals(i);    
    [ft_vals,om_vals]=FT_correlation(times,DO_func(times,om_1(omega,gamma),gamma).*resolution_func(times,sigma), n_points);
    dispname=sprintf('$n_{\\textrm{spline}} t = %.1f$',dt);
    plot(om_vals,abs(real(ft_vals)),...
        'DisplayName',str_points{i}, ...
        'Color',c_map(i,:), ...
        'LineWidth',2);
    hold on;
end
plot(om_vals,FFT_exact(om_vals,omega,gamma),...
        'DisplayName','exact', ...
        'Color',c_map(end,:), ...
        'LineWidth',2);
xlim([0,.2]);
ylim([1e-4,1e4]);
xlabel('$\omega$','Interpreter','latex');
ylabel('$|\hat{h}(\omega)|$','Interpreter','latex');
set(gca,'YScale','log')
legend show
set(legend,'interpreter','latex','location','northeast')

pos = get(gca, 'position');
dim = [pos(1)+.025, .1, pos(3), .3*pos(4)];
% str = {'$\omega_0 = 6 \cdot 10^{-2}$','$\sigma = 5\cdot 10^{3}$','$t_{\textrm{max}} = 10^4$'};
str = {omega_str,gamma_str,sigma_str,dt_str,tmax_str};
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'interpreter','latex');
% 
figname=sprintf('%s/DO_FFT_spline_add',fig_base);
fprintf('Completed figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end
