close all;
% clear all;
% load /data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/T_.15/samp_Dynamics.mat

fig_base = "/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/FFT_Testing/plots/FFT_Spline_freq_resolution";
saveswitch = 1;

i_q = 6; q_indices=i_q:numel(qbin):numel(averaging_times)*numel(qbin); 
q=qbin(i_q);
T=.15;
sqrtN=256;
om_step=.0005;
d_omega_vals=0:om_step:5*om_step;
c_map=linspecer(numel(d_omega_vals));
for i_om = 1:numel(d_omega_vals)
    d_omega = d_omega_vals(i_om);
    if (d_omega == 0)
        linestyle = '--';
        linewidth=1.5;
    else
        linestyle = '-';
        linewidth = 2.5;
%         linewidth = 2.1;
    end
    dispname=sprintf('$\\delta \\omega = %.1f \\times 10^{-3}$', 1e3*d_omega);
    [ft_vals{i_om},om_vals{i_om}]=FT_correlation(averaging_times,gmperpmperp(q_indices).*exp(-.5*d_omega^2.*averaging_times.^2), 0);
    plot(om_vals{i_om},real(ft_vals{i_om}),'DisplayName', ...
        dispname,'Color',c_map(i_om,:),...
        'LineWidth',linewidth,'LineStyle',linestyle);
    hold on; 
end
legend show;
set(legend,'interpreter','latex');
xlim([0,.05])
xlabel('$\omega$','interpreter','latex');
ylabel('$S_{m\perp}(q,\omega)$','interpreter','latex');

figname=sprintf('%s/resolution_yfull_T_%.3f_sqrtN_%d',fig_base,T,sqrtN);
fprintf('Completed figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

xlim([0,.1]);
ylim([0,1.3e4]);
figname=sprintf('%s/resolution_yzoom_T_%.3f_sqrtN_%d',fig_base,T,sqrtN);
fprintf('Completed figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end
