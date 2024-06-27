%% 0 Initialization
clear
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/TCF/NF_HighFreq',fig_base);

i_model = 1;
if (i_model == 1)
    curmodel="mxy";
    curtitle="MXY model";
    modelname="MXY";

%     sqrtN_vals = [16 32 64 128 256];
    sqrtN_vals = [16 32 64 128];
    dir_longtime="/data/scc/thobi/211201_LongerTime/mxy_3.00";
    dir_shorttime="/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00";

%     T_vals = [.11 .14 .17 .18];
%     T_dirs = {"T_.11" "T_.14" "T_.17" "T_.18"};
    
    T_vals = [.17];
    T_dirs = {"T_.17"};
%     T_vals = [.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52];
%     T_dirs = {"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
    
    sampfilename="samp_Dynamics";
    
    L_vals=[9.25,18.5,37,74,148];
    q_select=[1 2 3 4 5 6 7 8 10 11 13 14 15 18 21 22]; % Ignores multiple counts due to binning
    q_indices=[1 3 9 23];
    q_indices=[1 3 9];
    
    i_N_start = 2;
    i_N_max = 4;
        
elseif (i_model == 2)
        curmodel="xy";
        curtitle="XY model";

elseif (i_model == 3)
    curmodel="fmxy";
    curtitle="FMXY model";

    
elseif (i_model == 4)
    curmodel="xy_s";
    curtitle="XY S model";
    modelname="XY";

    sqrtN_vals = [16 32 64 128];
    dir_longtime="/data/scc/thobi/211201_LongerTime/xy_s";
    dir_shorttime="/data/scc/thobi/220201_ReducedSmapleStepDeltat/xy_s";

    T_vals = [.85 .91 .95 1.00];
    T_dirs = {"T_.85" "T_.91" "T_.95" "T_1.00"};
%     T_vals = [.17];
%     T_dirs = {"T_.17"};
%     T_vals = [.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52];
%     T_dirs = {"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
    
    sampfilename="samp_Dynamics";
    
    L_vals=sqrtN_vals;
%     q_select=[1 2 3 4 5 6 7 8 10 11 13 14 15 18 21 22]; % Ignores multiple counts due to binning
    q_indices=[1 3];
    i_N_start = 3;
    i_N_max = 4;
    
end
res_factor=30;
res_function = resolution_Gauss;
res_function=@(t,tau) resolution_Laplace_pleateau(t,tau,tau);

q_index_start = 6;

n_trunc = 4; % For NF Psi
z=1;

N_T = numel(T_vals);
N_N = numel(sqrtN_vals);
    

%% 1 Data collection

%% 2 The Plot
n_N = i_N_max - i_N_start + 1;
q_vals = zeros(1,n_N);
gamma_vals = zeros(1,n_N);
omega_1_vals = zeros(1,n_N);
om_max_vals = zeros(1,n_N);
n_period=4;
weightexp=1;
res_factor=[10,10,10,10];
res_factor=10*[1,1,1,1];
% res_factor=[8,5,3,3];
smooth_depth=5*[1,1,1,1];

% ylim_vec=[1e-3 Inf];
% eta=.24;
n_trunc = 5;

spline_factor=0;

xlim_vec=[1e-1 10*max((spline_factor)^(1/3),1)];

c_map=colormap_sqrtN();

eta_vals = zeros(1,N_T);
for i_T = 1:N_T
    T_dir=T_dirs{i_T};
    absM_vec=zeros(1,N_N);
    for i_N = 1:N_N
        sqrtN = sqrtN_vals(i_N);
        curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_longtime,sqrtN,T_dir,sampfilename);
        load(curfile,"absM_av");
        absM_vec(i_N) = absM_av;
    end
    eta_fitob = fit_eta_Magnetization_FS(absM_vec(1:N_N),L_vals(1:N_N));
    eta_vals(i_T) = eta_fitob.eta;
end

for i_T = 1:N_T
    figure
    T = T_vals(i_T);
    for i_N = i_N_start:i_N_max
        i = i_N - i_N_start + 1;
        i_q = q_indices(i);

        sqrtN = sqrtN_vals(i_N);
        L = L_vals(i_N);
        curdir_shorttime=sprintf('%s/sqrtN_%d/%s',dir_shorttime,sqrtN,T_dirs{i_T});
        curfile_shorttime=sprintf('%s/%s.mat',curdir_shorttime,sampfilename);
        curdir_longtime=sprintf('%s/sqrtN_%d/%s',dir_longtime,sqrtN,T_dirs{i_T});
        curfile_longtime=sprintf('%s/%s.mat',curdir_longtime,sampfilename);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SHORT TIME DATA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(curfile_shorttime,'averaging_times','gmperpmperp','qbin');

        q_vals(i) = qbin(i_q);
        t=averaging_times;
        cf=real(gmperpmperp(i_q:numel(qbin):end));

        
        c = fit_DampedOscillator_RealSpace(t,cf,n_period,weightexp,'omega_1');
        gamma_vals(i) = c(1);
        omega_1_vals(i) = c(2);


        tau=res_factor(i)/omega_1_vals(i) ;
%         tau=50/omega_1 ;
        res_vals=res_function(t,tau);
        
%         [ft_vals,om_vals]=FT_correlation(t, cf .* res_vals, 0);
        [ft_vals_shorttime,om_vals_shorttime]=FT_correlation(t, cf .* res_vals /cf(1), spline_factor * numel(cf));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LONG TIME DATA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(curfile_longtime,'averaging_times','gmperpmperp','qbin');

        q_vals(i) = qbin(i_q);
        t=averaging_times;
        cf=real(gmperpmperp(i_q:numel(qbin):end));

        
        c = fit_DampedOscillator_RealSpace(t,cf,n_period,weightexp,'omega_1');
        gamma_vals(i) = c(1);
        omega_1_vals(i) = c(2);


        tau=6/gamma_vals(i) ;
%         tau=50/omega_1 ;
        res_vals=res_function(t,tau);
%         res_vals = ones(size(t));
        
%         [ft_vals,om_vals]=FT_correlation(t, cf .* res_vals, 0);
        [ft_vals_longtime,om_vals_longtime]=FT_correlation(t, cf .* res_vals /cf(1), 0);


        [ft_max,i_max]=max(ft_vals_longtime);
        om_max_vals(i) = abs(om_vals_longtime(i_max));
        ft_max_vals(i) = ft_max;

%         ft_vals_supersmooth=smoothen(abs(real(ft_vals_shorttime)),1,smooth_depth(i));
%         ft_vals_minsmooth=smoothen(abs(real(ft_vals_shorttime)),1,1);
        non_smooth_ind = find( (om_vals_longtime > 0) .* (om_vals_longtime < 1.3 * om_max_vals(i)) );
        smooth_ind = find(om_vals_shorttime > 2 * om_max_vals(i));
        ft_vals_shorttime=smoothen(abs(real(ft_vals_shorttime(smooth_ind))),1,smooth_depth(i));
        om_vals_shorttime = om_vals_shorttime(smooth_ind);
        ft_vals_longtime=abs(real(ft_vals_longtime(non_smooth_ind)));
        om_vals_longtime = om_vals_longtime(non_smooth_ind);
%         combination_factor = ft_vals_longtime(non_smooth_ind(end) + 1) / ft_vals_supersmooth(smooth_ind(1));
%         combination_factor = 1;
%         ft_vals_combined = [ft_vals_longtime(non_smooth_ind),combination_factor*ft_vals_supersmooth(smooth_ind)];

%         om_vals_combined = [om_vals_longtime(non_smooth_ind),om_vals_shorttime(smooth_ind)];
%         ft_vals_shorttime = [ft_vals_minsmooth(non_smooth_ind),ft_vals_supersmooth(smooth_ind)];

%         dispname=sprintf('$N = (%d)^2$',sqrtN);
%         loglog(om_vals_combined / om_max_vals(i),om_max_vals(i)*real(ft_vals_combined)/cf(1),...
%             '-',...
%             'DisplayName',dispname,...
%             'LineWidth',1.5,...
%             'Color',c_map(i,:));
        dispname=sprintf('$N = (%d)^2$',sqrtN);
        hold on; 
        loglog(om_vals_longtime / om_max_vals(i),om_max_vals(i)*real(ft_vals_longtime),...
            '-',...
            'DisplayName',dispname,...
            'LineWidth',1.5,...
            'Color',c_map(i,:));        
        
        loglog(om_vals_shorttime / om_max_vals(i),om_max_vals(i)*real(ft_vals_shorttime),...
            '-',...
            'HandleVisibility','off',...
            'LineWidth',1.5,...
            'Color',c_map(i,:));

        loglog([om_vals_longtime(end),om_vals_shorttime(1)] / om_max_vals(i),om_max_vals(i)*[ft_vals_longtime(end),ft_vals_shorttime(1)],...
            ':',...
            'HandleVisibility','off',...
            'LineWidth',1.5,...
            'Color',c_map(i,:));
        


    end


    c_vals = omega_1_vals./q_vals;

    omega_1_cur = omega_1_vals(end);
    gamma_cur= gamma_vals(end);
    c_cur = c_vals(end);
    q_cur = q_vals(end);
    eta_cur = eta_vals(i_T);

    om_vals = logspace(log10(min(om_vals_longtime)),log10(max(om_vals_shorttime)),4e3);
    ft_vals=fitfunc_DO_reciprocal(om_vals,1,[gamma_cur,omega_1_cur]);
    S_q = sum(real(ft_vals(1:end))) * (om_vals_shorttime(2) - om_vals_shorttime(1));
    ft_scaled=c_cur * q_cur * real(ft_vals) / S_q;
    ft_scaled=c_cur * q_cur * real(ft_vals);
    loglog(om_vals / (c_cur * q_cur ), ft_scaled,...
            'Color','black', 'DisplayName','DO Fit', ... 'Color',c_map(i_N - i_N_start + 2,:),
            'LineStyle', '--', 'LineWidth',1.5, ...
            'Marker', 'none', 'MarkerSize', 5);
        
        
    y_vals=logspace(log10(xlim_vec(1)),log10(xlim_vec(2)),3e2);
%     y_vals=linspace(0,100,1e4);
    NF_Psi = NelsonFisher_Psi(y_vals,eta_cur,n_trunc);
    NF_Psi_scaled = NF_Psi; %c_cur * q_cur * 
    loglog(y_vals,NF_Psi_scaled,'-',...
            'Color','Black','DisplayName','NF', ...
            'LineWidth',1.5, ... 'Color',c_map(i_N - i_N_start + 3,:), 
            'Marker', 'none', 'MarkerSize', 5);



    xlim(xlim_vec); 
    ylim([.01*NF_Psi_scaled(end),3*max(ft_scaled)]);
    ylim([.01*NF_Psi_scaled(end),100]);
    hLegend=legend('Location','SouthWest','Interpreter','latex',...
    'NumColumns',1);
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'log','Yscale', 'log')
    
    hXLabel = xlabel('$\omega / cq$','interpreter','latex');
    hYLabel = ylabel('$cq S_{m\perp}(q,\omega)/\chi_{m\perp}(q)$','interpreter','latex');
    h_axis = gca;

    annotation_str = {modelname,sprintf('$N = (%d)^2$',sqrtN)};
    annotation_str = {modelname,sprintf('$n_{\\textrm{spline}} = %d$',spline_factor)};
    annotation_str = {modelname,sprintf('$T = %.3g$',T),sprintf('$q = %.3g$',q_cur)};
    dim=[.46 .16 .1 .2]; % Right of bottom legend
    dim=[.16 .43 .1 .2]; % Above left bottom legend
    h_annotation = annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
        'interpreter','latex',...
        'LineWidth', .5, ...
        'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
        'Color','black','FontSize', fontsize_annotation,...
        'BackgroundColor','white');

    
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    
    set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
    
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

    figname=sprintf('%s/%s_NFHighFreq_T_%.3g_q_%.3g_spline_factor_%d_smoothened',basedir,curmodel,T,q_cur,spline_factor);
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end

    
end


return
%% 
% S=load("/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00/sqrtN_128/T_.17/samp_Dynamics");
S=load("/data/scc/thobi/220201_ReducedSmapleStepDeltat/xy_s/sqrtN_128/T_.70/samp_Dynamics");

n_trunc=6;
eta_cur=-2*log(S.absM_av)/log(sqrt(2)*128);
y_vals=sort([logspace(-2,1,500),.9955:.001:1.01]);
NF_Psi=NelsonFisher_Psi(y_vals,eta_cur,n_trunc);

qbin=S.qbin;
gmperp=S.gmperpmperp;
t=S.averaging_times;
i_q=2^2*6;
q=qbin(i_q);
cf=real(gmperp(i_q:numel(qbin):end));

fitstr="exp(-gamma*x/2)*(cos(omega_1*x)+gamma/2/omega_1*sin(omega_1*x))";
fitfunc=@(t,ga,om) exp(-ga*t/2).*(cos(om*t)+ga/2/om*sin(om*t));
fitob=fit(t(:),cf(:)/cf(1),fitstr,'StartPoint',[.1,.1]);
om=fitob.omega_1;
ga=fitob.gamma;
c = om /q;

figure
plot(t,cf);
hold on;
plot(t,cf(1)*fitfunc(t,ga,om));

tau=1/ga;
res_vals=res_function(t,tau);
[ft_vals,om_vals]=FT_correlation(t, cf .* res_vals /cf(1), 1e6);

figure
ft_real=real(ft_vals);
ft_smooth=smoothen(abs(ft_real),1,3);
loglog(om_vals/om,ft_real);
hold on;
loglog(om_vals/om,ft_smooth);
loglog(y_vals,NF_Psi);
xlim([1e-1 4]);
