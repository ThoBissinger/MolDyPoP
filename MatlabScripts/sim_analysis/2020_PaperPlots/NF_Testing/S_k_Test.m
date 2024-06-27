%% 0 Initialization
run ../initialization_script
basedir=sprintf('%s/NF_Testing/S_k_Test',fig_base);
%% 1 Getting data
c=1e-2;
eta=.25;
r_vals=linspace(1e-1,148,1e3);
% t_vals=logspace(-2,6,1e3+1)';
t_vals=linspace(1e-1,1e5,1e4+1);
dr = r_vals(2) - r_vals(1);
dt = t_vals(2) - t_vals(1);
% r_vals=linspace(1e-2,1e2,5e2);
% t_vals=logspace(1e-2,1e3,1e3)';
r_mat=kron(r_vals,ones(size(t_vals)));
t_mat=kron(t_vals,ones(size(r_vals)));
y = c * t_mat ./ r_mat;
Phi = NelsonFisher_Phi(y,eta);
S_r = r_mat.^(-eta).*Phi;
[S_k,k_vals] = NelsonFisher_S_kt(r_vals,t_vals,c,eta);
dk = k_vals(2) - k_vals(1);
%% 2a Plot of S_k for fixed t, large r limit
figure
t_vals_plot=[0,3/c,10/c];
n_t = numel(t_vals_plot);
c_map = linspecer(n_t);
for i = 1:n_t
    t = t_vals_plot(i);
    i_t = find(t_vals >= t, 1);
    dispname = sprintf('$t = %d$', t);
    semilogx(k_vals,S_k(i_t,:), '-', ...
        'Color',c_map(i,:),'LineWidth',2,...
        'DisplayName',dispname);
    hold on;
    
end

% ylim([0, 1.2]);

h_axis = gca;
hXLabel = xlabel('$k$','interpreter','latex');
hYLabel = ylabel('$S(k,t)$','interpreter','latex');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legend, axes etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend = legend('Location', 'NorthEast','interpreter','latex','FontSize', fontsize_annotation,...
    'NumColumns',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Font
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)

annotation_str = {sprintf('$\\eta = %.2f$',eta), sprintf('$c = %d \\cdot 10^{-2}$',100*c)};
dim=[.43 .73 .1 .2];
annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

figname=sprintf('%s/Fixed_t_S_k',basedir);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');






%% 3a Plot of S_k for fixed k, focus on small t limit
figure
k_vals_plot=[dk,3*dk,5*dk];
n_k = numel(k_vals_plot);
c_map = linspecer(n_k);
for i = 1:n_k
    k = k_vals_plot(i);
    i_k = find(k_vals >= k, 1);
    dispname = sprintf('$k = %.2e$', k);
    semilogx(t_vals,real(S_k(:,i_k)), '-', ...
        'Color',c_map(i,:),'LineWidth',2,...
        'DisplayName',dispname);
    hold on;
end

% ylim([0, 1.2]);

h_axis = gca;
hXLabel = xlabel('$t$','interpreter','latex');
hYLabel = ylabel('$S(k,t)$','interpreter','latex');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legend, axes etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend = legend('Location', 'SouthWest','interpreter','latex','FontSize', fontsize_annotation,...
    'NumColumns',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Font
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)

annotation_str = {sprintf('$\\eta = %.2f$',eta), sprintf('$c = %d \\cdot 10^{-2}$',100*c)};
dim=[.15 .35 .1 .2];
annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

figname=sprintf('%s/Fixed_k_S_k',basedir);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');



%% 4 S_k_om test
figure
% om_vals_plot=[0,c/2,c,2*c];
% n_om = numel(om_vals_plot);
% c_map = linspecer(n_om);
% for i = 1:n_om
k_vals_plot=[5*dk,10*dk,20*dk];
n_k = numel(k_vals_plot);
c_map = linspecer(n_k);
for i = 1:n_k
    k = k_vals_plot(i);
    i_k = find(k_vals >= k, 1);
    [ft_vals,om_vals] = FT_correlation(t_vals,real(S_k(:,i_k)), 0);
    dispname = sprintf('$k = %.2e$', k);
    semilogx(om_vals / c / k,real(ft_vals), '-', ...
        'Color',c_map(i,:),'LineWidth',2,...
        'DisplayName',dispname);
    hold on;
end

% ylim([0, 1.2]);

h_axis = gca;
hXLabel = xlabel('$\omega/ck$','interpreter','latex');
hYLabel = ylabel('$S(k,t)$','interpreter','latex');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legend, axes etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend = legend('Location', 'SouthWest','interpreter','latex','FontSize', fontsize_annotation,...
    'NumColumns',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Font
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)

annotation_str = {sprintf('$\\eta = %.2f$',eta), sprintf('$c = %d \\cdot 10^{-2}$',100*c)};
dim=[.15 .35 .1 .2];
annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

figname=sprintf('%s/Fixed_k_S_k_om',basedir);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
