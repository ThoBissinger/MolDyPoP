%% 0 Initialization
run ../initialization_script
basedir=sprintf('%s/NF_Testing/S_k_om_Test',fig_base);
%% 1 Getting data
c=1e0;
eta=.25;
r_vals=linspace(1e-1,1e2,1e3);
t_vals=linspace(1e-1,1e4,3e3+1)';
dr = r_vals(2) - r_vals(1);
dt = t_vals(2) - t_vals(1);
% r_vals=linspace(1e-2,1e2,5e2);
% t_vals=logspace(1e-2,1e3,1e3)';
r_mat=kron(r_vals,ones(size(t_vals)));
t_mat=kron(t_vals,ones(size(r_vals)));
y = c * t_mat ./ r_mat;
Phi = NelsonFisher_Phi(y,eta);
S_r = r_mat.^(-eta).*Phi;
[S_k_om,k_vals,om_vals] = NelsonFisher_S_k_om(r_vals,t_vals,c,eta);
dk = k_vals(2) - k_vals(1);
dom = om_vals(2) - om_vals(1);
%% 2a Plot of S_k for fixed t, large r limit
figure
om_vals_plot=[0,c,3*c];
n_om = numel(om_vals_plot);
c_map = linspecer(n_om);
for i = 1:n_om
    om = om_vals_plot(i);
    i_om = find(om_vals >= om, 1);
    dispname = sprintf('$\\omega = %.3f$', om);
    semilogx(k_vals,real(S_k_om(i_om,:)), '-', ...
        'Color',c_map(i,:),'LineWidth',2,...
        'DisplayName',dispname);
    hold on;
    
end

% ylim([0, 1.2]);

h_axis = gca;
hXLabel = xlabel('$k$','interpreter','latex');
hYLabel = ylabel('$S(k,\omega)$','interpreter','latex');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legend, axes etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend = legend('Location', 'NorthEast','interpreter','latex','FontSize', fontsize_annotation,...
    'NumColumns',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Font
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)

annotation_str = {sprintf('$\\eta = %.2f$',eta), sprintf('$c = %d$',c)};
dim=[.43 .73 .1 .2];
annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');

figname=sprintf('%s/Fixed_om_S_k',basedir);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');






%% 3a Plot of S_k for fixed k, focus on small t limit
figure
k_vals_plot=[dk,10*dk,100*dk];
n_k = numel(k_vals_plot);
c_map = linspecer(n_k);
for i = 1:n_k
    k = k_vals_plot(i);
    i_k = find(k_vals >= k, 1);
    dispname = sprintf('$k = %.2e$', k);
    semilogx(om_vals,real(S_k_om(:,i_k)), '-', ...
        'Color',c_map(i,:),'LineWidth',2,...
        'DisplayName',dispname);
    hold on;
end

% ylim([0, 1.2]);

h_axis = gca;
hXLabel = xlabel('$\omega$','interpreter','latex');
hYLabel = ylabel('$S(k,\omega)$','interpreter','latex');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legend, axes etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend = legend('Location', 'NorthWest','interpreter','latex','FontSize', fontsize_annotation,...
    'NumColumns',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Font
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)

annotation_str = {sprintf('$\\eta = %.2f$',eta), sprintf('$c = %d$',c)};
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



