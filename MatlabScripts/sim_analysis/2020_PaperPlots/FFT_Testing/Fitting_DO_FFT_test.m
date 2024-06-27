% Fits damped oscillator data
close all;
% clear all;


if ispc
    fig_base = "C:/Users/tbiss/Thesis/MatlabScripts/sim_analysis/2020_PaperPlots/FFT_Testing/plots/Fitting_DO_FFT_test";
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\ExternalCodes\')
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\sim_analysis\FitFuncs\')
    addpath('C:\Users\tbiss\Thesis\MatlabScripts\sim_analysis\FourierTransform\')
else
    fig_base = "/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/FFT_Testing/plots/Fitting_DO_FFT_test";
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScriptsExternalCodes')
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScriptssim_analysis/FitFuncs')
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScriptssim_analysis/FourierTransform')
end
saveswitch = 1;

tmax=1e4;

fig_select=2;

DO_func=@(t,omega,gamma) exp(-.5*gamma * t) .* (cos(omega*t) + .5*gamma/omega*sin(omega*t));
FFT_exact=@(om,omega_0,gamma) 2*gamma*omega_0^2 ./ ...
    ( ( om.^2 - omega_0^2 ).^2 + om .^2 * gamma^2 );
resolution_gauss=@(t,sigma) exp(-.5 * t.^2 /sigma^2);
resolution_laplace=@(t,tau) exp(-abs(t)/tau);
om_1=@(omega_0,gamma) sqrt(omega_0^2 - .25*gamma^2);












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1 Fits to exact and FT data
if (ismember(1,fig_select))
    figure(1);

    dt=1e1;
    dt_str='$\Delta t = 10$';
    tmax=1e5;
    times=0:dt:tmax;
    tmax_str='$t_{\textrm{max}} = 10^5$';

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
        f_vals=DO_func(times,om_1(omega,gamma),gamma);
        [ft_vals_numeric,om_vals]=FT_correlation(times,f_vals, 0);
        ft_vals_exact=FFT_exact(om_vals,omega,gamma);
        fitob_numeric{i} = fit_DampedOscillator_Omega(om_vals,real(ft_vals_numeric),1);
        fitob_exact{i} = fit_DampedOscillator_Omega(om_vals,ft_vals_exact,1);
        plot(om_vals,fitob_numeric{i}.h_0 * FFT_exact(om_vals,fitob_numeric{i}.omega_0,fitob_numeric{i}.gamma) + fitob_numeric{i}.h_0 * fitob_numeric{i}.noiselevel,...
            'DisplayName',names{i}, ...
            'Color',c_map(i,:), ...
            'LineWidth',2);
        hold on;
        plot(om_vals,fitob_exact{i}.h_0 * FFT_exact(om_vals,fitob_exact{i}.omega_0,fitob_exact{i}.gamma) + fitob_exact{i}.h_0 * fitob_exact{i}.noiselevel,...
            'o', 'DisplayName',names{i}, ...
            'Color',c_map(i,:), ...
            'MarkerIndices',1:2:length(om_vals));
        fprintf("omega = %de-2, gamma = %de-4\nomega_num = %.4f, gamma_num = %.6f\nomega_ex = %.4f, gamma_ex = %.6f\n",omega*1e2, gamma*1e4,fitob_numeric{i}.omega_0,fitob_numeric{i}.gamma,fitob_exact{i}.omega_0,fitob_exact{i}.gamma)
    end
    legend show
    set(legend,'interpreter','latex','location','north')
    % set(gca,'YTick',[])
    xlabel('\omega');
    ylabel('$\hat{h}(\omega)$','Interpreter','latex');
    xlim([0,.05]);

    pos = get(gca, 'position');
    dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)];
    str = {dt_str,tmax_str};
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex');

    figname=sprintf('%s/Fit_compare_numeric_exact',fig_base);
    fprintf('Completed figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end

























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2 Plot of gamma_fit and omega_fit vs gamma
if (ismember(2,fig_select))
    figure(2);

    dt=1e1;
    dt_str='$\Delta t = 10$';
    tmax=1e5;
    times=0:dt:tmax;
    tmax_str='$t_{\textrm{max}} = 10^5$';

    omega=1e-3;
    om_str='$\omega_0 = 1 \cdot 10^{-3}$';
    n_gamma=30;
    gamma_vals=logspace(-5,-2,n_gamma);
    sigma=6e3;
    sig_str='$\sigma = 6 \cdot 10^{3}$';
    tau=2e4;
    tau_str='$\tau = 2 \cdot 10^{4}$';

    fitob_numeric = cell(1,numel(gamma_vals));
    fitob_res_gauss = cell(1,numel(gamma_vals));
    fitob_res_laplace = cell(1,numel(gamma_vals));
    fitob_exact = cell(1,numel(gamma_vals));
    
    gamma_numeric = zeros(1,numel(gamma_vals));
    gamma_res_gauss = zeros(1,numel(gamma_vals));
    gamma_res_laplace = zeros(1,numel(gamma_vals));
    gamma_exact = zeros(1,numel(gamma_vals));
    
    omega_numeric = zeros(1,numel(gamma_vals));
    omega_res_gauss = zeros(1,numel(gamma_vals));
    omega_res_laplace = zeros(1,numel(gamma_vals));
    omega_exact = zeros(1,numel(gamma_vals));
    for i = 1:numel(gamma_vals)
        gamma=gamma_vals(i);
    %     dispname=sprintf('$\\omega_0 = %d \\cdot 10^{-2},\\ \\gamma = %d \\cdot 10^{-4}$',omega*1e2,gamma*1e4);
        f_vals=DO_func(times,om_1(omega,gamma),gamma);
        [ft_vals_numeric,om_vals]=FT_correlation(times,f_vals, 0);
        [ft_vals_res_gauss,~]=FT_correlation(times,f_vals .* resolution_gauss(times,sigma), 0);
        [ft_vals_res_laplace,~]=FT_correlation(times,f_vals .* resolution_laplace(times,tau), 0);
        ft_vals_exact=FFT_exact(om_vals,omega,gamma);
        fitob_numeric{i} = fit_DampedOscillator_Omega(om_vals,real(ft_vals_numeric),1);
        fitob_res_gauss{i} = fit_DampedOscillator_Omega(om_vals,real(ft_vals_res_gauss),1);
        fitob_res_laplace{i} = fit_DampedOscillator_Omega(om_vals,real(ft_vals_res_laplace),1);
        fitob_exact{i} = fit_DampedOscillator_Omega(om_vals,ft_vals_exact,1);

        gamma_numeric(i)=fitob_numeric{i}.gamma;
        gamma_res_gauss(i)=fitob_res_gauss{i}.gamma;
        gamma_res_laplace(i)=fitob_res_laplace{i}.gamma;
        gamma_exact(i)=fitob_exact{i}.gamma;
        
        omega_numeric(i)=fitob_numeric{i}.omega_0;
        omega_res_gauss(i)=fitob_res_gauss{i}.omega_0;
        omega_res_laplace(i)=fitob_res_laplace{i}.omega_0;
        omega_exact(i)=fitob_exact{i}.omega_0;
    %     plot(om_vals,fitob_numeric{i}.h_0 * FFT_exact(om_vals,fitob_numeric{i}.omega_0,fitob_numeric{i}.gamma) + fitob_numeric{i}.h_0 * fitob_numeric{i}.noiselevel,...
    %         'DisplayName',names{i}, ...
    %         'Color',c_map(i,:), ...
    %         'LineWidth',2);
    %     hold on;
    %     plot(om_vals,fitob_exact{i}.h_0 * FFT_exact(om_vals,fitob_exact{i}.omega_0,fitob_exact{i}.gamma) + fitob_exact{i}.h_0 * fitob_exact{i}.noiselevel,...
    %         'o', 'DisplayName',names{i}, ...
    %         'Color',c_map(i,:), ...
    %         'MarkerIndices',1:2:length(om_vals));
    %     fprintf("omega = %de-2, gamma = %de-4\nomega_num = %.4f, gamma_num = %.6f\nomega_ex = %.4f, gamma_ex = %.6f\n",omega*1e2, gamma*1e4,fitob_numeric{i}.omega_0,fitob_numeric{i}.gamma,fitob_exact{i}.omega_0,fitob_exact{i}.gamma)
    end

    c_map=linspecer(4);
    loglog(gamma_vals,gamma_exact,'o',...
        'Displayname','direct fit',...
        'Color',c_map(1,:), 'LineWidth',2);
    hold on;
    loglog(gamma_vals,gamma_numeric,'o',...
        'Displayname','numeric, no resolution',...
        'Color',c_map(2,:), 'LineWidth',2);
    loglog(gamma_vals,gamma_res_gauss,'o',...
        'Displayname','Gauss resolution',...
        'Color',c_map(3,:), 'LineWidth',2);
    loglog(gamma_vals,gamma_res_laplace,'o',...
        'Displayname','Laplace resolution',...
        'Color',c_map(4,:), 'LineWidth',2);
    loglog(gamma_vals,gamma_vals,'--',...
        'Displayname','$\gamma$ exact',...
        'Color','black', 'LineWidth',2);
    legend show
    set(legend,'interpreter','latex','location','northwest')
    % set(gca,'YTick',[])
    xlabel('\gamma');
    ylabel('$\gamma_{\textrm{fit}}$','Interpreter','latex');
    xlim([gamma_vals(1),gamma_vals(end)]);

    pos = get(gca, 'position');
    dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)]; %top right
    dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
    dim = [pos(3)-.05, .1, pos(3), .3*pos(4)]; % bottom right
    str = {dt_str,tmax_str,om_str,sig_str,tau_str};
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex');

    figname=sprintf('%s/Fit_compare_gamma_gamma',fig_base);
    fprintf('Completed figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    
    figure(102)
    c_map=linspecer(4);
    loglog(gamma_vals,omega_exact,'o',...
        'Displayname','direct fit',...
        'Color',c_map(1,:), 'LineWidth',2);
    hold on;
    loglog(gamma_vals,omega_numeric,'o',...
        'Displayname','numeric, no resolution',...
        'Color',c_map(2,:), 'LineWidth',2);
    loglog(gamma_vals,omega_res_gauss,'o',...
        'Displayname','Gauss resolution',...
        'Color',c_map(3,:), 'LineWidth',2);
    loglog(gamma_vals,omega_res_laplace,'o',...
        'Displayname','Laplace resolution',...
        'Color',c_map(4,:), 'LineWidth',2);
%     loglog(gamma_vals,omega*ones(size(gamma_vals)),'--',...
%         'Displayname','$\omega$ exact',...
%         'Color','black', 'LineWidth',2);
    legend show
    set(legend,'interpreter','latex','location','northwest')
    % set(gca,'YTick',[])
    xlabel('\gamma');
    ylabel('$\gamma_{\textrm{fit}}$','Interpreter','latex');
    xlim([gamma_vals(1),gamma_vals(end)]);

    pos = get(gca, 'position');
    dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)]; %top right
    dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
    dim = [pos(3)-.05, .1, pos(3), .3*pos(4)]; % bottom right
    str = {dt_str,tmax_str,om_str,sig_str,tau_str};
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex');

    figname=sprintf('%s/Fit_compare_gamma_omega',fig_base);
    fprintf('Completed figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3 Plot of omega_fit and gamma_fit vs omega
if (ismember(3,fig_select))
    figure(3);

    dt=1e1;
    dt_str='$\Delta t = 10$';
    tmax=1e5;
    times=0:dt:tmax;
    tmax_str='$t_{\textrm{max}} = 10^5$';

    omega_vals=logspace(-5,-1,30);
    n_omega=30;
    gamma=5e-4;
    om_str='$\gamma = 5 \cdot 10^{-4}$';
    sigma=6e3;
    sig_str='$\sigma = 6 \cdot 10^{3}$';
    tau=2e4;
    tau_str='$\tau = 2 \cdot 10^{4}$';
    
    fitob_numeric = cell(1,numel(omega_vals));
    fitob_res_gauss = cell(1,numel(omega_vals));
    fitob_res_laplace = cell(1,numel(omega_vals));
    fitob_exact = cell(1,numel(omega_vals));
    
    gamma_numeric = zeros(1,numel(omega_vals));
    gamma_res_gauss = zeros(1,numel(omega_vals));
    gamma_res_laplace = zeros(1,numel(omega_vals));
    gamma_exact = zeros(1,numel(omega_vals));

    omega_numeric = zeros(1,numel(omega_vals));
    omega_res_gauss = zeros(1,numel(omega_vals));
    omega_res_laplace = zeros(1,numel(omega_vals));
    omega_exact = zeros(1,numel(omega_vals));
    for i = 1:numel(omega_vals)
        omega=omega_vals(i);
    %     dispname=sprintf('$\\omega_0 = %d \\cdot 10^{-2},\\ \\gamma = %d \\cdot 10^{-4}$',omega*1e2,gamma*1e4);
        f_vals=DO_func(times,om_1(omega,gamma),gamma);
        [ft_vals_numeric,om_vals]=FT_correlation(times,f_vals, 0);
        [ft_vals_res_gauss,~]=FT_correlation(times,f_vals .* resolution_gauss(times,sigma), 0);
        [ft_vals_res_laplace,~]=FT_correlation(times,f_vals .* resolution_laplace(times,tau), 0);
        ft_vals_exact=FFT_exact(om_vals,omega,gamma);
        fitob_numeric{i} = fit_DampedOscillator_Omega(om_vals,real(ft_vals_numeric),1);
        fitob_res_gauss{i} = fit_DampedOscillator_Omega(om_vals,real(ft_vals_res_gauss),1);
        fitob_res_laplace{i} = fit_DampedOscillator_Omega(om_vals,real(ft_vals_res_laplace),1);
        fitob_exact{i} = fit_DampedOscillator_Omega(om_vals,ft_vals_exact,1);

        gamma_numeric(i)=fitob_numeric{i}.gamma;
        gamma_res_gauss(i)=fitob_res_gauss{i}.gamma;
        gamma_res_laplace(i)=fitob_res_laplace{i}.gamma;
        gamma_exact(i)=fitob_exact{i}.gamma;
        
        omega_numeric(i)=fitob_numeric{i}.omega_0;
        omega_res_gauss(i)=fitob_res_gauss{i}.omega_0;
        omega_res_laplace(i)=fitob_res_laplace{i}.omega_0;
        omega_exact(i)=fitob_exact{i}.omega_0;
    %     plot(om_vals,fitob_numeric{i}.h_0 * FFT_exact(om_vals,fitob_numeric{i}.omega_0,fitob_numeric{i}.gamma) + fitob_numeric{i}.h_0 * fitob_numeric{i}.noiselevel,...
    %         'DisplayName',names{i}, ...
    %         'Color',c_map(i,:), ...
    %         'LineWidth',2);
    %     hold on;
    %     plot(om_vals,fitob_exact{i}.h_0 * FFT_exact(om_vals,fitob_exact{i}.omega_0,fitob_exact{i}.gamma) + fitob_exact{i}.h_0 * fitob_exact{i}.noiselevel,...
    %         'o', 'DisplayName',names{i}, ...
    %         'Color',c_map(i,:), ...
    %         'MarkerIndices',1:2:length(om_vals));
    %     fprintf("omega = %de-2, gamma = %de-4\nomega_num = %.4f, gamma_num = %.6f\nomega_ex = %.4f, gamma_ex = %.6f\n",omega*1e2, gamma*1e4,fitob_numeric{i}.omega_0,fitob_numeric{i}.gamma,fitob_exact{i}.omega_0,fitob_exact{i}.gamma)
    end

    c_map=linspecer(4);
    loglog(omega_vals,omega_exact,'o',...
        'Displayname','direct fit',...
        'Color',c_map(1,:), 'LineWidth',2);
    hold on;
    loglog(omega_vals,omega_numeric,'o',...
        'Displayname','numeric, no resolution',...
        'Color',c_map(2,:), 'LineWidth',2);
    loglog(omega_vals,omega_res_gauss,'o',...
        'Displayname','Gauss resolution',...
        'Color',c_map(3,:), 'LineWidth',2);
    loglog(omega_vals,omega_res_laplace,'o',...
        'Displayname','Laplace resolution',...
        'Color',c_map(4,:), 'LineWidth',2);
    loglog(omega_vals,omega_vals,'--',...
        'Displayname','$\gamma_0$ exact',...
        'Color','black', 'LineWidth',2);
    legend show
    set(legend,'interpreter','latex','location','northwest')
    % set(gca,'YTick',[])
    xlabel('\omega_0');
    ylabel('$\omega_{0,\textrm{fit}}$','Interpreter','latex');
    xlim([omega_vals(1),omega_vals(end)]);

    pos = get(gca, 'position');
    dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)]; %top right
    dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
    dim = [pos(3)-.05, .15, pos(3), .3*pos(4)]; % bottom right
    str = {dt_str,tmax_str,om_str,sig_str,tau_str};
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex');

    figname=sprintf('%s/Fit_compare_omega_omega',fig_base);
    fprintf('Completed figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    % gammafit vs omega
    figure(103);
    c_map=linspecer(4);
    loglog(omega_vals,gamma_exact,'o',...
        'Displayname','direct fit',...
        'Color',c_map(1,:), 'LineWidth',2);
    hold on;
    loglog(omega_vals,gamma_numeric,'o',...
        'Displayname','numeric, no resolution',...
        'Color',c_map(2,:), 'LineWidth',2);
    loglog(omega_vals,gamma_res_gauss,'o',...
        'Displayname','Gauss resolution',...
        'Color',c_map(3,:), 'LineWidth',2);
    loglog(omega_vals,gamma_res_laplace,'o',...
        'Displayname','Laplace resolution',...
        'Color',c_map(4,:), 'LineWidth',2);
%     loglog(omega_vals,gamma*ones(size(omega_vals)),'--',...
%         'Displayname','$\gamma_0$ exact',...
%         'Color','black', 'LineWidth',2);
    legend show
    set(legend,'interpreter','latex','location','northwest')
    % set(gca,'YTick',[])
    xlabel('\omega_0');
    ylabel('$\gamma_{\textrm{fit}}$','Interpreter','latex');
    xlim([omega_vals(1),omega_vals(end)]);

    pos = get(gca, 'position');
    dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)]; %top right
    dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
    dim = [pos(3)-.05, .15, pos(3), .3*pos(4)]; % bottom right
    str = {dt_str,tmax_str,om_str,sig_str,tau_str};
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex');

    figname=sprintf('%s/Fit_compare_omega_gamma',fig_base);
    fprintf('Completed figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4 Plot of gamma_fit vs gamma for different Gaussian resolution
if (ismember(4,fig_select))
    figure(4);

    
    dt=1e1;
    dt_str='$\Delta t = 10$';
    tmax=1e5;
    times=0:dt:tmax;
    tmax_str='$t_{\textrm{max}} = 10^5$';

    n_gamma=15;
    gamma_vals=logspace(-5,-2,n_gamma);
    n_sigma=5;
    omega=2e-2;
    om_str='$\omega_0 = 2 \cdot 10^{-2}$';

    sigma_vals=logspace(3,5,n_sigma);

    fitob=cell(n_sigma,n_gamma);
    gamma_vec = zeros(n_sigma,n_gamma);
    for i = 1:n_gamma
        gamma=gamma_vals(i);
    %     dispname=sprintf('$\\omega_0 = %d \\cdot 10^{-2},\\ \\gamma = %d \\cdot 10^{-4}$',omega*1e2,gamma*1e4);
        f_vals=DO_func(times,om_1(omega,gamma),gamma);
        for j = 1:n_sigma
            sigma=sigma_vals(j);
            [ft_vals,om_vals]=FT_correlation(times,f_vals .* resolution_gauss(times,sigma), 0);
            fitob{j,i} = fit_DampedOscillator_Omega(om_vals,real(ft_vals),1);
            gamma_vec(j,i)=fitob{j,i}.gamma;
        end
    
    end

    c_map=linspecer(n_sigma+1);
    for j=1:n_sigma
        sigma=sigma_vals(j);
        dispname=sprintf('$\\sigma = %.1e $',sigma);
        loglog(gamma_vals,gamma_vec(j,:),'o',...
            'Displayname',dispname,...
            'Color',c_map(j,:), 'LineWidth',2);
        hold on;
    end
    loglog(gamma_vals,gamma_vals,'--',...
        'Displayname','$\gamma$ exact',...
        'Color','black', 'LineWidth',2);
    legend show
    set(legend,'interpreter','latex','location','northwest')
    % set(gca,'YTick',[])
    xlabel('\gamma');
    ylabel('$\gamma_{\textrm{fit}}$','Interpreter','latex');
    xlim([gamma_vals(1),gamma_vals(end)]);
    ylim([1e-5,1e-1])

    pos = get(gca, 'position');
    dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)]; %top right
    dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
    dim = [pos(3)-.05, .15, pos(3), .3*pos(4)]; % bottom right
    str = {dt_str,tmax_str,om_str};
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex');

    figname=sprintf('%s/Fit_resolution_compare_gauss_gamma',fig_base);
    fprintf('Completed figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5 Plot of omega_fit vs omega for different Gaussian resolution
if (ismember(5,fig_select))
    figure(5);

    dt=1e1;
    dt_str='$\Delta t = 10$';
    tmax=1e5;
    times=0:dt:tmax;
    tmax_str='$t_{\textrm{max}} = 10^5$';

    n_omega=15;
    n_sigma=5;
    omega_vals=logspace(-5,-1,n_omega);
    sigma_vals=logspace(3,5,n_sigma);
    gamma=5e-4;
    gamma_str='$\gamma = 5 \cdot 10^{-4}$';
%     sigma=6e3;
%     sig_str='$\sigma = 6 \cdot 10^{3}$';
%     tau=2e4;
%     tau_str='$\tau = 2 \cdot 10^{4}$';

    fitob=cell(n_sigma,n_omega);
    omega_vec = zeros(n_sigma,n_omega);
    for i = 1:n_omega
        omega=omega_vals(i);
    %     dispname=sprintf('$\\omega_0 = %d \\cdot 10^{-2},\\ \\gamma = %d \\cdot 10^{-4}$',omega*1e2,gamma*1e4);
        f_vals=DO_func(times,om_1(omega,gamma),gamma);
        for j = 1:n_sigma
            sigma=sigma_vals(j);
            [ft_vals,om_vals]=FT_correlation(times,f_vals .* resolution_gauss(times,sigma), 0);
            fitob{j,i} = fit_DampedOscillator_Omega(om_vals,real(ft_vals),1);
            omega_vec(j,i)=fitob{j,i}.omega_0;
        end
    
    end

    c_map=linspecer(n_sigma+1);
    for j=1:n_sigma
        sigma=sigma_vals(j);
        dispname=sprintf('$\\sigma = %.1e $',sigma);
        loglog(omega_vals,omega_vec(j,:),'o',...
            'Displayname',dispname,...
            'Color',c_map(j,:), 'LineWidth',2);
        hold on;
    end
    loglog(omega_vals,omega_vals,'--',...
        'Displayname','$\omega_0$ exact',...
        'Color','black', 'LineWidth',2);
    legend show
    set(legend,'interpreter','latex','location','northwest')
    % set(gca,'YTick',[])
    xlabel('\omega_0');
    ylabel('$\omega_{0,\textrm{fit}}$','Interpreter','latex');
    xlim([omega_vals(1),omega_vals(end)]);

    pos = get(gca, 'position');
    dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)]; %top right
    dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
    dim = [pos(3)-.05, .15, pos(3), .3*pos(4)]; % bottom right
    str = {dt_str,tmax_str,gamma_str};
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex');

    figname=sprintf('%s/Fit_resolution_compare_gauss_omega',fig_base);
    fprintf('Completed figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6 Plot of gamma_fit vs gamma for different Laplacian resolution
if (ismember(6,fig_select))
    figure(6);

    
    dt=1e1;
    dt_str='$\Delta t = 10$';
    tmax=1e5;
    times=0:dt:tmax;
    tmax_str='$t_{\textrm{max}} = 10^5$';

    n_gamma=15;
    gamma_vals=logspace(-5,-2,n_gamma);
    n_tau=5;
    omega=2e-2;
    om_str='$\omega_0 = 2 \cdot 10^{-2}$';

    tau_vals=logspace(3,5,n_tau);

    fitob=cell(n_tau,n_gamma);
    gamma_vec = zeros(n_tau,n_gamma);
    for i = 1:n_gamma
        gamma=gamma_vals(i);
    %     dispname=sprintf('$\\omega_0 = %d \\cdot 10^{-2},\\ \\gamma = %d \\cdot 10^{-4}$',omega*1e2,gamma*1e4);
        f_vals=DO_func(times,om_1(omega,gamma),gamma);
        for j = 1:n_tau
            tau=tau_vals(j);
            [ft_vals,om_vals]=FT_correlation(times,f_vals .* resolution_laplace(times,tau), 0);
            fitob{j,i} = fit_DampedOscillator_Omega(om_vals,real(ft_vals),1);
            gamma_vec(j,i)=fitob{j,i}.gamma;
        end
    
    end

    c_map=linspecer(n_tau+1);
    for j=1:n_tau
        tau=tau_vals(j);
        dispname=sprintf('$\\tau = %.1e $',tau);
        loglog(gamma_vals,gamma_vec(j,:),'o',...
            'Displayname',dispname,...
            'Color',c_map(j,:), 'LineWidth',2);
        hold on;
    end
    loglog(gamma_vals,gamma_vals,'--',...
        'Displayname','$\gamma$ exact',...
        'Color','black', 'LineWidth',2);
    legend show
    set(legend,'interpreter','latex','location','northwest')
    % set(gca,'YTick',[])
    xlabel('\gamma');
    ylabel('$\gamma_{\textrm{fit}}$','Interpreter','latex');
    xlim([gamma_vals(1),gamma_vals(end)]);
    ylim([1e-5,1e-1])

    pos = get(gca, 'position');
    dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)]; %top right
    dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
    dim = [pos(3)-.05, .15, pos(3), .3*pos(4)]; % bottom right
    str = {dt_str,tmax_str,om_str};
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex');

    figname=sprintf('%s/Fit_resolution_compare_laplace_gamma',fig_base);
    fprintf('Completed figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7 Plot of omega_fit vs omega for different Laplacian resolution
if (ismember(7,fig_select))
    figure(7);

    dt=1e1;
    dt_str='$\Delta t = 10$';
    tmax=1e5;
    times=0:dt:tmax;
    tmax_str='$t_{\textrm{max}} = 10^5$';

    n_omega=15;
    n_tau=5;
    omega_vals=logspace(-5,-1,n_omega);
    tau_vals=logspace(3,5,n_tau);
    gamma=5e-4;
    gamma_str='$\gamma = 5 \cdot 10^{-4}$';
%     tau=6e3;
%     sig_str='$\tau = 6 \cdot 10^{3}$';
%     tau=2e4;
%     tau_str='$\tau = 2 \cdot 10^{4}$';

    fitob=cell(n_tau,n_omega);
    omega_vec = zeros(n_tau,n_omega);
    for i = 1:n_omega
        omega=omega_vals(i);
    %     dispname=sprintf('$\\omega_0 = %d \\cdot 10^{-2},\\ \\gamma = %d \\cdot 10^{-4}$',omega*1e2,gamma*1e4);
        f_vals=DO_func(times,om_1(omega,gamma),gamma);
        for j = 1:n_tau
            tau=tau_vals(j);
            [ft_vals,om_vals]=FT_correlation(times,f_vals .* resolution_laplace(times,tau), 0);
            fitob{j,i} = fit_DampedOscillator_Omega(om_vals,real(ft_vals),1);
            omega_vec(j,i)=fitob{j,i}.omega_0;
        end
    
    end

    c_map=linspecer(n_tau+1);
    for j=1:n_tau
        tau=tau_vals(j);
        dispname=sprintf('$\\tau = %.1e $',tau);
        loglog(omega_vals,omega_vec(j,:),'o',...
            'Displayname',dispname,...
            'Color',c_map(j,:), 'LineWidth',2);
        hold on;
    end
    loglog(omega_vals,omega_vals,'--',...
        'Displayname','$\omega_0$ exact',...
        'Color','black', 'LineWidth',2);
    legend show
    set(legend,'interpreter','latex','location','northwest')
    % set(gca,'YTick',[])
    xlabel('\omega_0');
    ylabel('$\omega_{0,\textrm{fit}}$','Interpreter','latex');
    xlim([omega_vals(1),omega_vals(end)]);

    pos = get(gca, 'position');
    dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)]; %top right
    dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
    dim = [pos(3)-.05, .15, pos(3), .3*pos(4)]; % bottom right
    str = {dt_str,tmax_str,gamma_str};
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex');

    figname=sprintf('%s/Fit_resolution_compare_laplace_omega',fig_base);
    fprintf('Completed figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8 Plot of gamma_fit vs gamma for different noise levels
if (ismember(8,fig_select))
    figure(8);
    
    dt=1e1;
    dt_str='$\Delta t = 10$';
    tmax=1e5;
    times=0:dt:tmax;
    tmax_str='$t_{\textrm{max}} = 10^5$';

    n_gamma=15;
    gamma_vals=logspace(-4,-2,n_gamma);
    n_noise=5;
    omega=2e-2;
    om_str='$\omega_0 = 2 \cdot 10^{-2}$';

    noise_vals=logspace(-2,0,n_noise);
%     markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
    markers = {'+','o','x','s','*','p','h'};

    fitob=cell(n_noise,n_gamma);
    gamma_vec = zeros(n_noise,n_gamma);
    for i = 1:n_gamma
        gamma=gamma_vals(i);
    %     dispname=sprintf('$\\omega_0 = %d \\cdot 10^{-2},\\ \\gamma = %d \\cdot 10^{-4}$',omega*1e2,gamma*1e4);
        for j = 1:n_noise
            noiselevel=noise_vals(j);
            f_vals=DO_func(times,om_1(omega,gamma),gamma) + noiselevel*(rand(size(times))-.5);
            [ft_vals,om_vals]=FT_correlation(times,f_vals, 0);
            fitob{j,i} = fit_DampedOscillator_Omega(om_vals,real(ft_vals),1);
            gamma_vec(j,i)=fitob{j,i}.gamma;
        end
    
    end

    c_map=linspecer(n_noise+1);
    for j=1:n_noise
        noiselevel=noise_vals(j);
        dispname=sprintf('$\\xi = %.1e $',noiselevel);
        loglog(gamma_vals,gamma_vec(j,:),markers{j},...
            'Displayname',dispname,...
            'Color',c_map(j,:), 'LineWidth',1.2);
        hold on;
    end
    loglog(gamma_vals,gamma_vals,'--',...
        'Displayname','$\gamma$ exact',...
        'Color','black', 'LineWidth',2);
    legend show
    set(legend,'interpreter','latex','location','northwest')
    % set(gca,'YTick',[])
    xlabel('\gamma');
    ylabel('$\gamma_{\textrm{fit}}$','Interpreter','latex');
    xlim([gamma_vals(1),gamma_vals(end)]);
    ylim([1e-5,1e-1])

    pos = get(gca, 'position');
    dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)]; %top right
    dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
    dim = [pos(3)-.05, .15, pos(3), .3*pos(4)]; % bottom right
    str = {dt_str,tmax_str,om_str};
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex');

    figname=sprintf('%s/Fit_noise_compare_gamma',fig_base);
    fprintf('Completed figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9 Plot of omega_fit vs omega for different noise levels
if (ismember(9,fig_select))
    figure(9);

    dt=1e1;
    dt_str='$\Delta t = 10$';
    tmax=1e5;
    times=0:dt:tmax;
    tmax_str='$t_{\textrm{max}} = 10^5$';

    n_omega=15;
    n_noise=5;
    omega_vals=logspace(-4,-1,n_omega);
    noise_vals=logspace(-3,-1,n_noise);
%     markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
    markers = {'+','o','x','s','*','p','h'};
    gamma=5e-4;
    gamma_str='$\gamma = 5 \cdot 10^{-4}$';
%     tau=6e3;
%     sig_str='$\tau = 6 \cdot 10^{3}$';
%     tau=2e4;
%     tau_str='$\tau = 2 \cdot 10^{4}$';

    fitob=cell(n_noise,n_omega);
    omega_vec = zeros(n_noise,n_omega);
    for i = 1:n_omega
        omega=omega_vals(i);
    %     dispname=sprintf('$\\omega_0 = %d \\cdot 10^{-2},\\ \\gamma = %d \\cdot 10^{-4}$',omega*1e2,gamma*1e4);
        for j = 1:n_noise
            noiselevel=noise_vals(j);
            f_vals=DO_func(times,om_1(omega,gamma),gamma) + noiselevel*(rand(size(times))-.5);
            [ft_vals,om_vals]=FT_correlation(times,f_vals, 0);
            fitob{j,i} = fit_DampedOscillator_Omega(om_vals,real(ft_vals),1);
            omega_vec(j,i)=fitob{j,i}.omega_0;
        end
    
    end

    c_map=linspecer(n_noise+1);
    for j=1:n_noise
        noiselevel=noise_vals(j);
        dispname=sprintf('$\\xi = %.1e $',noiselevel);
        loglog(omega_vals,omega_vec(j,:),markers{j},...
            'Displayname',dispname,...
            'Color',c_map(j,:), 'LineWidth',1.2);
        hold on;
    end
    loglog(omega_vals,omega_vals,'--',...
        'Displayname','$\omega_0$ exact',...
        'Color','black', 'LineWidth',2);
    legend show
    set(legend,'interpreter','latex','location','northwest')
    % set(gca,'YTick',[])
    xlabel('\omega_0');
    ylabel('$\omega_{0,\textrm{fit}}$','Interpreter','latex');
    xlim([omega_vals(1),omega_vals(end)]);

    pos = get(gca, 'position');
    dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)]; %top right
    dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
    dim = [pos(3)-.05, .15, pos(3), .3*pos(4)]; % bottom right
    str = {dt_str,tmax_str,gamma_str};
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'interpreter','latex');

    figname=sprintf('%s/Fit_noise_compare_omega',fig_base);
    fprintf('Completed figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end

