%% 0 Initialization
run ../initialization_script

%% 1
n_t = 1e3+1;
n_spline_vals = [4:6];
t_max=1e4;

t=linspace(0,t_max,n_t);
omega_0=4e-2;
gamma=1e-3;
omega_1 = om_1(omega_0,gamma);

figure
tau=20/gamma;
resolution_vals=resolution_Laplace(t,tau);
resolution_vals=resolution_Gauss(t,3e3);
% resolution_vals=ones(size(t));

cf=fitfunc_DO(t,1,[gamma,omega_1]);
c_map=linspecer(numel(n_spline_vals));
for i_spline = 1:numel(n_spline_vals)
    n_spline = n_spline_vals(i_spline);

    [ft_vals,om_vals]=FT_correlation(t, cf .* resolution_vals, 10^n_spline);
    S_q = sum(real(ft_vals(1:end))) * (om_vals(2) - om_vals(1));
    
    dispname = sprintf('$ n_{\\textrm{spline}} = 10^{%d}$',n_spline);
    loglog(om_vals,real(ft_vals) ,...
        'Color',c_map(i_spline,:),'DisplayName',dispname, ...
        'LineStyle', '-', 'LineWidth',2, ...
        'Marker', '.', 'MarkerSize', 5);
    hold on;
end
xlim([0,Inf]);
ylim([1e-2/max(om_vals), 2/gamma])
ft_exact=fitfunc_DO_reciprocal(om_vals,1,[gamma,omega_0]);
S_q = sum(abs(ft_exact(1:end))) * (om_vals(2) - om_vals(1));
loglog(om_vals,abs(ft_exact) ,...
    'Color','Black','DisplayName','$(\omega L^z)^{-1}$', ...
    'LineStyle', '--', 'LineWidth',2, ...
    'Marker', 'none', 'MarkerSize', 5);

h_axis = gca;
hXLabel = xlabel('$\omega L^z$','interpreter','latex');
hYLabel = ylabel('$S_{m}(q,\omega) / L^z S_{m}(q)$','interpreter','latex');
hLegend = legend('Location', 'Northeast','interpreter','latex','FontSize', fontsize_annotation,...
    'NumColumns',1);
set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)