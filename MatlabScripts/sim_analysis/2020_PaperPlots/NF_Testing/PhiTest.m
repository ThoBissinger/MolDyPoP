%% 0 Initialization
run ../initialization_script
basedir=sprintf('%s/NF_Testing/PhiTest',fig_base);
%% 1 Getting data
% r_vals=linspace(1e-1,148,5e2);
% t_vals=linspace(0,1e3,1e3+1);
r_vals=logspace(-2,3,5e2);
t_vals=logspace(-2,5,1e3)';
r_mat=kron(r_vals,ones(size(t_vals)));
t_mat=kron(t_vals,ones(size(r_vals)));
c=1e-2;
eta=.25;
y = c * t_mat ./ r_mat;
Phi = NelsonFisher_Phi(y,eta);

%% 2a Plot of Phi for fixed t, large r limit
figure
t_vals_plot=[1e1,1e2,5e2];
n_t = numel(t_vals_plot);
c_map = linspecer(n_t);
for i = 1:n_t
    t = t_vals_plot(i);
    i_t = find(t_vals >= t, 1);
    dispname = sprintf('$t = %d$', t);
    semilogx(r_vals,Phi(i_t,:), '-', ...
        'Color',c_map(i,:),'LineWidth',2,...
        'DisplayName',dispname);
    hold on;
end

ylim([0, 1.2]);

h_axis = gca;
hXLabel = xlabel('$r$','interpreter','latex');
hYLabel = ylabel('$r^{\eta} S(r,t)$','interpreter','latex');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legend, axes etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend = legend('Location', 'SouthEast','interpreter','latex','FontSize', fontsize_annotation,...
    'NumColumns',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Font
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)

annotation_str = {sprintf('$\\eta = %.2f$',eta), sprintf('$c = %d \\cdot 10^{-2}$',100*c)};
dim=[.435 .138 .1 .2];
annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

figname=sprintf('%s/Fixed_t_phi_large_r',basedir);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');



%% 2b Plot of Phi for fixed t, small r limit
figure
t_vals_plot=[1e1,1e2,5e2];
n_t = numel(t_vals_plot);
c_map = linspecer(n_t);
for i = 1:n_t
    t = t_vals_plot(i);
    i_t = find(t_vals >= t, 1);
    dispname = sprintf('$t = %d$', t);
    semilogx(r_vals,(2 * c * t)^eta * Phi(i_t,:)./r_vals.^eta, '-', ...
        'Color',c_map(i,:),'LineWidth',2,...
        'DisplayName',dispname);
    hold on;
end
h_axis = gca;
hXLabel = xlabel('$r$','interpreter','latex');
hYLabel = ylabel('$2 ct S(r,t)$','interpreter','latex');




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
dim=[.37 .138 .1 .2];
annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

figname=sprintf('%s/Fixed_t_phi_small_r',basedir);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');





%% 3a Plot of Phi for fixed r, focus on small t limit
figure
r_vals_plot=[1e0,5e0,1e1];
n_r = numel(t_vals_plot);
c_map = linspecer(n_r);
for i = 1:n_r
    r = r_vals_plot(i);
    i_r = find(r_vals >= r, 1);
    dispname = sprintf('$r = %d$', r);
    semilogx(t_vals,Phi(:,i_r), '-', ...
        'Color',c_map(i,:),'LineWidth',2,...
        'DisplayName',dispname);
    hold on;
end

ylim([0, 1.2]);

h_axis = gca;
hXLabel = xlabel('$t$','interpreter','latex');
hYLabel = ylabel('$r^{\eta} S(r,t)$','interpreter','latex');



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
dim=[.375 .138 .1 .2];
annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

figname=sprintf('%s/Fixed_r_phi_small_t',basedir);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');






%% 3b Plot of Phi for fixed r, focus on large t limit
figure
r_vals_plot=[1e0,2e0,5e0];
n_r = numel(t_vals_plot);
c_map = linspecer(n_r);
for i = 1:n_r
    r = r_vals_plot(i);
    i_r = find(r_vals >= r, 1);
    dispname = sprintf('$r = %d$', r);
    semilogx(t_vals,(2 * c * t_vals).^eta .* Phi(:,i_r)/r^eta, '-', ...
        'Color',c_map(i,:),'LineWidth',2,...
        'DisplayName',dispname);
    hold on;
end

ylim([0, 1.2]);

h_axis = gca;
hXLabel = xlabel('$t$','interpreter','latex');
hYLabel = ylabel('$2 ct S(r,t)$','interpreter','latex');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legend, axes etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend = legend('Location', 'NorthWest','interpreter','latex','FontSize', fontsize_annotation,...
    'NumColumns',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Font
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)

annotation_str = {sprintf('$\\eta = %.2f$',eta), sprintf('$c = %d \\cdot 10^{-2}$',100*c)};
dim=[.15 .5 .1 .2];
annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

figname=sprintf('%s/Fixed_r_phi_large_t',basedir);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');

