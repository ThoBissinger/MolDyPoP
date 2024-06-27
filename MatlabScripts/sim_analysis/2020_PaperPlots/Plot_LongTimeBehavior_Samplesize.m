%% 0 Loading
initialization_script
matfile_collect='/data/scc/thobi/211201_LongerTime/fmxy/sqrtN_16/T_.17/samp_Dynamics_collect.mat';
matfile_compact='/data/scc/thobi/211201_LongerTime/fmxy/sqrtN_16/T_.17/samp_Dynamics.mat';
load(matfile_collect,'gmperpmperp_collect');
load(matfile_compact,'averaging_times','qbin');




runmax=3e3;
T=.17;
sqrtN=16;
n_q=numel(qbin);
n_t=numel(averaging_times);
i_q=1;
q=qbin(i_q);
q_indices=i_q:n_q:n_q*n_t;


%% 1 Plot
n_run=numel(run_num_vec);
run_num_vec=[47, 188, 750 3e3];
c_map = linspecer(n_run);
t_min=700;
t_max=1000;
figure
for i_run = 1:n_run
    run_num=run_num_vec(i_run);
    run_indices=runmax-run_num+1:runmax;
    cf=mean(gmperpmperp_collect(run_indices,q_indices));
    cf=real(cf)/real(cf(1));
    dispname=sprintf('$n_{\\textrm{samp}} = %d$',run_num);
    plot(averaging_times,cf,'-',...
        'Color',c_map(i_run,:),...
        'LineWidth',2,...
        'DisplayName',dispname);
    hold on;
    yline(1/sqrt(run_num),'Color',c_map(i_run,:),...
        'LineWidth',2.3,...
        'HandleVisibility','off');
end
xlim([t_min,t_max]);
ymax=1/sqrt(min(run_num_vec));
ylim([-ymax,2.5*ymax]);

h_axis = gca;
hXLabel = xlabel('$t$','interpreter','latex');
hYLabel = ylabel('$C_{m\perp}(q,t)/\chi_{m\perp}(q)$','interpreter','latex');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legend, axes etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegend = legend('Location', 'northeast','interpreter','latex','FontSize', fontsize_annotation,...
    'NumColumns',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Font
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)

annotation_str = {sprintf('$N = (%d)^2$',sqrtN),sprintf('$T = %.3f$',T),...
    sprintf('$q = %.3f$',q)};
dim=[.16 .7 .1 .2];
annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
    'interpreter','latex',...
    'VerticalAlignment','middle', 'HorizontalAlignment','left',...
    'Color','black','FontSize', fontsize_annotation,...
    'BackgroundColor','white');
figname='plots/LongTimeBehavior_Samplesize/sqrtN_16_T_.17';
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');



%% Plot average fluctuation at large t vs n_samp
samp_step=10;
n_step=runmax/samp_step;
samp_vec=samp_step:samp_step:runmax;
var_vec=zeros(1,n_step);
absmean_vec=zeros(1,n_step);
i_start=find(averaging_times > t_min, 1);
for i = 1:n_step
    runmin = runmax-samp_step*i + 1;
    run_indices=runmin:runmax;
    cf=mean(gmperpmperp_collect(run_indices,q_indices));
    cf=real(cf(i_start:end))/real(cf(1));
    
    var_vec(i) = var(cf);
    absmean_vec(i) = mean(abs(cf));
end
figure
plot(samp_vec,var_vec);
hold on;
plot(samp_vec,absmean_vec);
