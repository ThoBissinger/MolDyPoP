%% 1 Different system sizes and different temperatures
clear;
close all;
initialization_script;

fontsize_annotation = 10;
fontsize_axis = 8;
fontsize_ax_labels = 10;
fontsize_labels = 14;
basedir=sprintf('%s/plots/Magnetization_Evolution/ang_MSD',fig_base);

% dataset = "AdjustedTime";
dataset = "LinearTime";

t_scaling = "t";
% t_scaling = "tbyL";
% t_scaling = "tbyL2";

y_scaling = "lin";
y_scaling = "Lto2-eta";

if dataset == "AdjustedTime"
    sqrtN_vals = [16, 32, 64, 128];
    T_vals = [.11, .14, .17, .185];
    T_dirs = {"T_.11", "T_.14", "T_.17", "T_.185"};
    dirs=sprintfc("/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_%d/",sqrtN_vals);
elseif dataset == "LinearTime"
    sqrtN_vals = [16, 32, 64, 128, 256];
    T_vals = [.14, .17, .19, .21];
    T_dirs = {"T_.14", "T_.17", "T_.19", "T_.20"};
    dirs=sprintfc("/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_%d/",sqrtN_vals);
end
dirs_logtime=sprintfc("/data/scc/thobi/210727_LepriRuffo_GoodResolution/mxy_3.00/sqrtN_%d/",sqrtN_vals);    
% T_dirs=sprintfc("T_%.3g",T_vals); 

if t_scaling == "t"
    t_scale=@(t,L) t;
    x_str = "$t$";
elseif t_scaling == "tbyL"
    t_scale=@(t,L) t/L;
    x_str = "$t/L$";
elseif t_scaling == "tbyL2"
    t_scale=@(t,L) t/L^2;
    x_str = "$t/L^2$";
end



L_vals = [9.25, 18.5, 37, 74, 148];
fig_indices = [1 2 4 5];
sampfilename = "samp_Dynamics_M_TimeEvolution";
sampfilename_logtime = "samp_LepriRuffo_M_TimeEvolution";
N_T = numel(T_vals);
N_N = numel(sqrtN_vals);
figure
for i_T = 1:N_T
    T = T_vals(i_T);
    c_map=linspecer(N_N);
    for i_N = 1:N_N
        sqrtN = sqrtN_vals(i_N);
        L = L_vals(i_N);
        curdir = sprintf("%s/%s",dirs{i_N},T_dirs{i_T});
        curdir_logtime = sprintf("%s/%s",dirs_logtime{i_N},T_dirs{i_T});
        load(sprintf('%s/%s',curdir,sampfilename),...
            "averaging_times","M_ang_MSD");%,"Mx_collect","My_collect");
        t_lintime=averaging_times;
        y_lintime=M_ang_MSD;
        load(sprintf('%s/%s',curdir_logtime,sampfilename_logtime),...
            "averaging_times","M_ang_MSD");
        t_logtime=averaging_times;
        y_logtime=M_ang_MSD;

        [t,y]=combine_cf(t_logtime,y_logtime,t_lintime,y_lintime);
%         absM_av = mean(mean(sqrt(Mx_collect.^2 + My_collect.^2)));
        t_max = max(averaging_times);
        n_t = numel(averaging_times);
        t_max_exp = round(log10(t_max));
        t_max_noexp = t_max / 10^t_max_exp;
%         subfig_index = (i_N - 1) *N_T + i_T;
%         subplot(N_N,N_T,subfig_index)
        subfig_index = fig_indices(i_T);
        subplot(2,3,subfig_index)
        dispname=sprintf('$N = (%d)^2$',sqrtN);
        loglog(t_scale(t,L),y,...
            'Color',c_map(i_N,:),...
            'Displayname',dispname,...
            'LineWidth',2);
        hold on;

    
    
        
        set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'XGrid', 'off', 'YGrid', 'off', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
            'XTick',logspace(-10,4,15),'YTick',logspace(-10,4,15),...
            'LineWidth', .5,...
            'XScale','log','YScale','log')

%         if (i_N == 1)
%             h_Tlabel = text(.5, 1.1,sprintf("$T = %.3f$",T)','units','normalized',...
%                 'interpreter','latex',...
%                 'HorizontalAlignment','center',...
%                 'fontsize',fontsize_labels);
%         end
        
    end 
    h_axis = gca;
    h_title = title(sprintf("$T = %.3f$",T),'interpreter','latex');
    hXLabel = xlabel(x_str,'interpreter','latex');
    hYLabel = ylabel('$\langle \Delta \phi_m^2(t)\rangle$','interpreter','latex');
    
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    
    set(h_title,'FontName','cmr12','fontsize',fontsize_labels)
    xlim([max(t_scale(5e-1,L),1e-3) Inf]);
    ylim([1e-5,Inf])
%     xx=logspace(0,1);
%     plot(t_scale(xx,L),1e-6*xx.^2,'--',...
%         'LineWidth',2,...
%         'Color','black');
%     xx=logspace(2.5,4);
%     plot(t_scale(xx,L),1e-6*xx,'--',...
%         'LineWidth',2,...
%         'Color','black');
end

subplot (2, 3, [3 6]) % merge remaining subplots and put legend here
for i_N = 1:N_N
    dispname=sprintf('$N = (%d)^2$', sqrtN_vals(i_N));
    plot(nan, nan,'Color',c_map(i_N,:),'DisplayName',dispname,...
        'LineWidth',2) % plot nans (hack to generate correct legend but plot no data)
    hold on;
end
% hLegend = legend(sprintfc('$N = (%d)^2$', sqrtN_vals), 'Location', 'west',...
%     'interpreter','latex'); 
axis off
hLegend = legend('interpreter','latex','Location','West');
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

set(gcf,'units','centimeters','OuterPosition',[0 0 pagewidth_cm .5*pageheight_cm]);


figname=sprintf('%s/mxy_M_ang_MSD_%s_%s',basedir,dataset,t_scaling);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

%% 2 Different temperatures
clear;
close all;
initialization_script;

fontsize_annotation = 10;
fontsize_axis = 8;
fontsize_ax_labels = 10;
fontsize_labels = 14;
basedir=sprintf('%s/plots/Magnetization_Evolution/ang_MSD',fig_base);

dataset = "AdjustedTime";
% dataset = "LinearTime";

t_scaling = "t";
% t_scaling = "tbyL";
% t_scaling = "tbyL2";
if dataset == "AdjustedTime"
    sqrtN_vals = [16 32 64 128];
    L_vals = [9.25 18.5 37 74];
%     T_vals = [.11, .14, .165, .167, .169, .17, .171, .173, .175, .18, .185];
%     T_dirs = {"T_.11", "T_.14", "T_.165", "T_.167", "T_.169", "T_.17", ...
%         "T_.171", "T_.173", "T_.175", "T_.18", "T_.185"};
    T_vals = [.11, .14, .165, .17, .175, .18, .185];
    T_dirs = {"T_.11", "T_.14", "T_.165", "T_.17", ...
        "T_.175", "T_.18", "T_.185"};
    dirs=sprintfc("/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_%d/",sqrtN_vals);
    dirs_logtime=sprintfc("/data/scc/thobi/210727_LepriRuffo_GoodResolution/mxy_3.00/sqrtN_%d/",sqrtN_vals);    

%     t_longtime_fct=@(sqrtN) 1e4*(128/sqrtN)^2;
elseif dataset == "LinearTime"
    sqrtN_vals = 128;
    L_vals = [74];
    T_vals = [.14, .17, .19, .21];
    dirs=sprintfc("/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_%d/",sqrtN_vals);
end
% T_dirs=sprintfc("T_%.3g",T_vals); 

if t_scaling == "t"
    t_scale=@(t,L) t;
    x_str = "$t$";
elseif t_scaling == "tbyL"
    t_scale=@(t,L) t/L;
    x_str = "$t/L$";
elseif t_scaling == "tbyL2"
    t_scale=@(t,L) t/L^2;
    x_str = "$t/L^2$";
end
sampfilename = "samp_Dynamics_M_TimeEvolution";
sampfilename_logtime = "samp_LepriRuffo_M_TimeEvolution";
% sampfilename = "samp_Dynamics";
N_T = numel(T_vals);
N_N = numel(sqrtN_vals);

fitob_shorttime_cell = cell(N_N,N_T);
fitob_longtime_cell = cell(N_N,N_T);
a_shorttime_vec = zeros(N_N,N_T);
b_shorttime_vec = zeros(N_N,N_T);
D_longtime_vec = zeros(N_N,N_T);
b_longtime_vec = zeros(N_N,N_T);
c_longtime_vec = zeros(N_N,N_T);

fitfunc_shorttime =@(a,b,x) a*x.^b;
% fitfunc_longtime =@(D,b,c,x) D*x.^b + c;
fitfunc_longtime =@(D,b,x) D*x.^b;

for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    figure
    c_map=turbo(N_T+1);
    c_map=c_map(2:end,:);
%     c_map=cool(N_T);
%     N_T_small = sum(T_vals <= .17);
%     c_map=[cool(N_T_small);summer(N_T-N_T_small)];
%     c_map=othercolor('RdYlGn11',N_T);
%     c_map=othercolor('Set13',N_T);

    for i_T = 1:N_T
        T = T_vals(i_T);
        L = L_vals(i_N);
        curdir = sprintf("%s/%s",dirs{i_N},T_dirs{i_T});
        curdir_logtime = sprintf("%s/%s",dirs_logtime{i_N},T_dirs{i_T});
        load(sprintf('%s/%s',curdir,sampfilename),...
            "averaging_times","M_ang_MSD");%,"Mx_collect","My_collect");
        t_lintime=averaging_times;
        y_lintime=M_ang_MSD;
        load(sprintf('%s/%s',curdir_logtime,sampfilename_logtime),...
            "averaging_times","M_ang_MSD");
        t_logtime=averaging_times;
        y_logtime=M_ang_MSD;

        [t,y]=combine_cf(t_logtime,y_logtime,t_lintime,y_lintime);

%         absM_av = mean(mean(sqrt(Mx_collect.^2 + My_collect.^2)));
        t_max = max(averaging_times);
        n_t = numel(averaging_times);
        t_max_exp = round(log10(t_max));
        t_max_noexp = t_max / 10^t_max_exp;
        dispname=sprintf('$T = %.3g$',T);
        loglog(t_scale(t,L),y,'-',...
            'Color',c_map(i_T,:),...
            'MarkerSize',3,...
            'Displayname',dispname,...
            'LineWidth',2);
        hold on;

    
    
        
        set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'XGrid', 'off', 'YGrid', 'off', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
            'XTick',logspace(-10,4,15),'YTick',logspace(-10,4,15),...
            'LineWidth', .5,...
            'XScale','log','YScale','log')

        i_max_shorttime = find(t > .2, 1);
        t_shorttime=t(1:i_max_shorttime);
        y_shorttime=y(1:i_max_shorttime);
        init_b = 2;
        init_a = y_shorttime(end)/t_shorttime(end)^init_b;
    
        fitob_shorttime = fit(t_shorttime(:),y_shorttime(:),fittype(fitfunc_shorttime),...
            'StartPoint',[init_a,init_b]);
%         [init_a, init_b; fitob_shorttime.a, fitob_shorttime.b]
       
        
        i_min_longtime = find(t > t(end)/5,1);
        t_longtime=t(i_min_longtime:end);
        y_longtime=y(i_min_longtime:end);
        init_D = (y_longtime(end) - y_longtime(1))/...
            (t_longtime(end) - t_longtime(1));
        init_b = 1;
        init_c = y_longtime(1) - init_D * t_longtime(1);
%         fitob_longtime = fit(t_longtime(:),y_longtime(:),fittype(fitfunc_longtime),...
%             'StartPoint',[init_D,init_b,init_c]);
        fitob_longtime = fit(t_longtime(:),y_longtime(:),fittype(fitfunc_longtime),...
            'StartPoint',[init_D,init_b]);
    %     fitob_longtime = fit(t_longtime(:),y_longtime(:),fitfunc_longtime,...
    
        a_shorttime_vec(i_N,i_T) = fitob_shorttime.a;
        b_shorttime_vec(i_N,i_T) = fitob_shorttime.b;
        
        D_longtime_vec(i_N,i_T) = fitob_longtime.D;
        b_longtime_vec(i_N,i_T) = fitob_longtime.b;
%         c_longtime_vec(i_N,i_T) = fitob_longtime.c;
        
        fitob_shorttime_cell{i_N,i_T} = fitob_shorttime;
        fitob_longtime_cell{i_N,i_T} = fitob_longtime;

%         tt = t_shorttime;
%         a_cur = a_shorttime_vec(i_N,i_T);
%         b_cur = b_shorttime_vec(i_N,i_T);
%         y_fit = fitfunc_shorttime(a_cur,b_cur,tt);
%         plot(tt,y_fit,'--',...
%             'Color','black','LineWidth',2,...
%             'HandleVisibility','off')
%     
%         tt = t_longtime;
%         D_cur = D_longtime_vec(i_N,i_T);
%         b_cur = b_longtime_vec(i_N,i_T);
% %         c_cur = c_longtime_vec(i_N,i_T);
%         y_fit = fitfunc_longtime(D_cur,b_cur,tt);
% %         y_fit = fitfunc_longtime(D_cur,b_cur,c_cur,tt);
%         plot(tt,y_fit,'--',...
%             'Color','black','LineWidth',2,...
%             'HandleVisibility','off')
%         
    end 
    xlim([1e-1, Inf]);
    h_axis = gca;
%     h_title = title(sprintf("$T = %.3f$",T),'interpreter','latex');
    hXLabel = xlabel(x_str,'interpreter','latex');
    hYLabel = ylabel('$\langle \Delta \phi_m^2(t)\rangle$','interpreter','latex');
    
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    
    tt = logspace(-1,0);
    a_cur = 2*a_shorttime_vec(i_N,end);
    b_cur = b_shorttime_vec(i_N,end);
    y_fit = fitfunc_shorttime(a_cur,b_cur,tt);
    plot(tt,y_fit,'--',...
        'Color','black','LineWidth',2,...
        'HandleVisibility','off')
    annotation_str = sprintf('$\\propto t^{%.2f}$',b_cur);
    h_tshort_annotation = text(tt(end), y_fit(end),annotation_str,...
        'interpreter','latex',...
        'HorizontalAlignment','right',...
        'VerticalAlignment','bottom',...
        'fontsize',fontsize_labels,...
        'Color','black');

    tt = logspace(3,4);
    D_cur = D_longtime_vec(i_N,3);
    b_cur = b_longtime_vec(i_N,3);
%     c_cur = c_longtime_vec(i_N,1);
    y_fit = fitfunc_longtime(D_cur,b_cur,tt);
    plot(tt,y_fit,'--',...
        'Color','black','LineWidth',2,...
        'HandleVisibility','off')
    annotation_str = sprintf('$\\propto t^{%.2f}$',b_cur);
    h_tlong_annotation = text(tt(end/2), y_fit(end/2),annotation_str,...
        'interpreter','latex',...
        'HorizontalAlignment','left',...
        'VerticalAlignment','top',...
        'fontsize',fontsize_labels,...
        'Color','black');
%     set(h_title,'FontName','cmr12','fontsize',fontsize_labels)
%     xlim([t_scale(t(1),L) Inf]);
%     xx=logspace(0,1);
%     plot(t_scale(xx,L),1e-6*xx.^2,'--',...
%         'LineWidth',2,...
%         'Color','black');
%     xx=logspace(2.5,4);
%     plot(t_scale(xx,L),1e-6*xx,'--',...
%         'LineWidth',2,...
%         'Color','black');

    hLegend = legend('interpreter','latex','Location','SouthEast',...
        'NumColumns',3);
    set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
    
    set(gcf,'units','centimeters','OuterPosition',[0 0 pagewidth_cm .5*pageheight_cm]);
    
    
    figname=sprintf('%s/mxy_M_ang_MSD_sqrtN_%d_%s_%s',basedir,sqrtN,dataset,t_scaling);
    
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end

end

%% 2.1 Short-time coefficient vs L
figure
hold on;
c_map=turbo(N_T+1);
c_map=c_map(2:end,:);
fitfunc_L =@(a0,bL,x) a0*x.^bL;
fitfunc_L_log =@(a0,bL,x) log(a0)+bL*log(x);
fit_a0_vec = zeros(size(T_vals));
fit_bL_a_vec = zeros(size(T_vals));
for i_T = 1:N_T
    T = T_vals(i_T);
    dispname=sprintf('$T = %.3g$',T);
    loglog(L_vals,a_shorttime_vec(:,i_T),'s--',...
        'Color',c_map(i_T,:),...
        'MarkerSize',6,...
        'Displayname',dispname,...
        'LineWidth',2);

    init_bL_a = -2;
    init_a0 = a_shorttime_vec(1,i_T) / L_vals(1)^init_bL_a;
    fitob_L = fit(L_vals(:),log(a_shorttime_vec(:,i_T)),fittype(fitfunc_L_log),...
        'StartPoint',[init_a0,init_bL_a],...
        'Lower',[0 -Inf]);
%     fitob_L = fit(L_vals(:),a_shorttime_vec(:,i_T),fittype(fitfunc_L),...
%         'StartPoint',[init_a0,init_bL]);
    fit_a0_vec(i_T) = fitob_L.a0;
    fit_bL_a_vec(i_T) = fitob_L.bL;
    annotation_str = sprintf('$\\ \\propto L^{%.2f}$',fit_bL_a_vec(i_T));
    h_tshort_annotation = text(L_vals(end), a_shorttime_vec(end,i_T),annotation_str,...
        'interpreter','latex',...
        'HorizontalAlignment','left',...
        'VerticalAlignment','baseline',...
        'fontsize',fontsize_labels,...
        'Color','black');

end
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'off', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'XTick',logspace(-10,4,15),'YTick',logspace(-10,4,15),...
    'LineWidth', .5,...
    'XScale','log','YScale','log')
h_axis = gca;
hXLabel = xlabel('L','interpreter','latex');
hYLabel = ylabel('short time $a$','interpreter','latex');

set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);

hLegend = legend('interpreter','latex','Location','NorthEast',...
        'NumColumns',3);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

set(gcf,'units','centimeters','OuterPosition',[0 0 pagewidth_cm .5*pageheight_cm]);


figname=sprintf('%s/mxy_M_ang_MSD_a_coeff',basedir);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end




















%% 2.2 Long-time coefficient vs L
figure
hold on;
c_map=turbo(N_T+1);
c_map=c_map(2:end,:);
fitfunc_L =@(a0,bL,x) a0*x.^bL;
fitfunc_L_log =@(a0,bL,x) log(a0)+bL*log(x);
fit_D0_vec = zeros(size(T_vals));
fit_bL_D_vec = zeros(size(T_vals));
for i_T = 1:N_T
    T = T_vals(i_T);
    dispname=sprintf('$T = %.3g$',T);
    loglog(L_vals,D_longtime_vec(:,i_T),'s--',...
        'Color',c_map(i_T,:),...
        'MarkerSize',6,...
        'Displayname',dispname,...
        'LineWidth',2);

    init_bL_D = -1;
    init_D0 = D_longtime_vec(1,i_T) / L_vals(1)^init_bL_D;
    fitob_L = fit(L_vals(:),log(D_longtime_vec(:,i_T)),fittype(fitfunc_L_log),...
        'StartPoint',[init_D0,init_bL_D],...
        'Lower',[0 -Inf]);
%     fitob_L = fit(L_vals(:),a_shorttime_vec(:,i_T),fittype(fitfunc_L),...
%         'StartPoint',[init_a0,init_bL]);
    fit_D0_vec(i_T) = fitob_L.a0;
    fit_bL_D_vec(i_T) = fitob_L.bL;
    annotation_str = sprintf('$\\ \\propto L^{%.2f}$',fit_bL_D_vec(i_T));
    h_tshort_annotation = text(L_vals(end), D_longtime_vec(end,i_T),annotation_str,...
        'interpreter','latex',...
        'HorizontalAlignment','left',...
        'VerticalAlignment','baseline',...
        'fontsize',fontsize_labels,...
        'Color','black');

end
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XGrid', 'off', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'XTick',logspace(-10,4,15),'YTick',logspace(-10,4,15),...
    'LineWidth', .5,...
    'XScale','log','YScale','log')
h_axis = gca;
hXLabel = xlabel('L','interpreter','latex');
hYLabel = ylabel('long time $D$','interpreter','latex');
set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);



hLegend = legend('interpreter','latex','Location','NorthEast',...
        'NumColumns',3);
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

set(gcf,'units','centimeters','OuterPosition',[0 0 pagewidth_cm .5*pageheight_cm]);


figname=sprintf('%s/mxy_M_ang_MSD_D_coeff',basedir);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end



%% 3 Magnetization Decay
clear;
close all;
initialization_script;

fontsize_annotation = 10;
fontsize_axis = 8;
fontsize_ax_labels = 10;
fontsize_labels = 14;
basedir=sprintf('%s/plots/Magnetization_Evolution/ang_MSD',fig_base);

dataset = "AdjustedTime";
% dataset = "LinearTime";

t_scaling = "t";
% t_scaling = "tbyL";
% t_scaling = "tbyL2";
if dataset == "AdjustedTime"
    sqrtN_vals = 128;
    L_vals = [74];
%     T_vals = [.11, .14, .165, .167, .169, .17, .171, .173, .175, .18, .185];
%     T_dirs = {"T_.11", "T_.14", "T_.165", "T_.167", "T_.169", "T_.17", ...
%         "T_.171", "T_.173", "T_.175", "T_.18", "T_.185"};
    T_vals = [.11, .14, .165, .17, .175, .18, .185];
    T_dirs = {"T_.11", "T_.14", "T_.165", "T_.17", ...
        "T_.175", "T_.18", "T_.185"};
    dirs=sprintfc("/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_%d/",sqrtN_vals);
    dirs_logtime=sprintfc("/data/scc/thobi/210727_LepriRuffo_GoodResolution/mxy_3.00/sqrtN_%d/",sqrtN_vals);    

%     t_longtime_fct=@(sqrtN) 1e4*(128/sqrtN)^2;
elseif dataset == "LinearTime"
    sqrtN_vals = 128;
    L_vals = [74];
    T_vals = [.14, .17, .19, .21];
    dirs=sprintfc("/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_%d/",sqrtN_vals);
end
% T_dirs=sprintfc("T_%.3g",T_vals); 

if t_scaling == "t"
    t_scale=@(t,L) t;
    x_str = "$t$";
elseif t_scaling == "tbyL"
    t_scale=@(t,L) t/L;
    x_str = "$t/L$";
elseif t_scaling == "tbyL2"
    t_scale=@(t,L) t/L^2;
    x_str = "$t/L^2$";
end
sampfilename = "samp_Dynamics_M_TimeEvolution";
sampfilename_logtime = "samp_LepriRuffo_M_TimeEvolution";
% sampfilename = "samp_Dynamics";
N_T = numel(T_vals);
N_N = numel(sqrtN_vals);

fitob_shorttime_cell = cell(N_N,N_T);
fitob_longtime_cell = cell(N_N,N_T);
a_shorttime_vec = zeros(N_N,N_T);
b_shorttime_vec = zeros(N_N,N_T);
D_longtime_vec = zeros(N_N,N_T);
b_longtime_vec = zeros(N_N,N_T);
c_longtime_vec = zeros(N_N,N_T);

fitfunc_shorttime =@(a,b,x) a*x.^b;
% fitfunc_longtime =@(D,b,c,x) D*x.^b + c;
fitfunc_longtime =@(D,b,x) D*x.^b;

for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    figure
    c_map=turbo(N_T+1);
    c_map=c_map(2:end,:);
%     c_map=cool(N_T);
%     N_T_small = sum(T_vals <= .17);
%     c_map=[cool(N_T_small);summer(N_T-N_T_small)];
%     c_map=othercolor('RdYlGn11',N_T);
%     c_map=othercolor('Set13',N_T);

    for i_T = 1:N_T
        T = T_vals(i_T);
        L = L_vals(i_N);
        curdir = sprintf("%s/%s",dirs{i_N},T_dirs{i_T});
        curdir_logtime = sprintf("%s/%s",dirs_logtime{i_N},T_dirs{i_T});
        load(sprintf('%s/%s',curdir,sampfilename),...
            "averaging_times","M_ang_MSD","Mx_collect","My_collect");
        t_lintime=averaging_times;
        y_lintime=mean(Mx_collect(:,1).*Mx_collect + My_collect(:,1).*My_collect);

        load(sprintf('%s/%s',curdir_logtime,sampfilename_logtime),...
            "averaging_times","M_ang_MSD","Mx_collect","My_collect");
        t_logtime=averaging_times;
        y_logtime=mean(Mx_collect(:,1).*Mx_collect + My_collect(:,1).*My_collect);

        [t,y]=combine_cf(t_logtime,y_logtime,t_lintime,y_lintime);

%         absM_av = mean(mean(sqrt(Mx_collect.^2 + My_collect.^2)));
        t_max = max(averaging_times);
        n_t = numel(averaging_times);
        t_max_exp = round(log10(t_max));
        t_max_noexp = t_max / 10^t_max_exp;
        dispname=sprintf('$T = %.3g$',T);
        loglog(t_scale(t,L),y,'-',...
            'Color',c_map(i_T,:),...
            'MarkerSize',3,...
            'Displayname',dispname,...
            'LineWidth',2);
        hold on;

    
    
        
        set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'XGrid', 'off', 'YGrid', 'off', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
            'XTick',logspace(-10,4,15),'YTick',0:.1:1,...
            'LineWidth', .5,...
            'XScale','log','YScale','linear')

    end 
    ylim([0 .61]);
    h_axis = gca;
%     h_title = title(sprintf("$T = %.3f$",T),'interpreter','latex');
    hXLabel = xlabel(x_str,'interpreter','latex');
    hYLabel = ylabel('$\langle \textbf{m}\cdot\textbf{m}(t)\rangle$','interpreter','latex');
    
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    

    hLegend = legend('interpreter','latex','Location','SouthEast',...
        'NumColumns',3);
    set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)
    
    set(gcf,'units','centimeters','OuterPosition',[0 0 pagewidth_cm .5*pageheight_cm]);
    
    
    figname=sprintf('%s/mxy_M_SCF_sqrtN_%d_%s_%s',basedir,sqrtN,dataset,t_scaling);
    
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end

end




%% 4 Different system sizes, different temperatures, y-axis scaled with L^{2-eta}
clear;
close all;
initialization_script;

fontsize_annotation = 10;
fontsize_axis = 8;
fontsize_ax_labels = 10;
fontsize_labels = 14;
basedir=sprintf('%s/plots/Magnetization_Evolution/ang_MSD',fig_base);

dataset = "AdjustedTime";
% dataset = "LinearTime";

t_scaling = "t";
% t_scaling = "tbyL";
% t_scaling = "tbyL2";

y_scaling = "lin";
y_scaling = "Lto2-eta";

if dataset == "AdjustedTime"
    sqrtN_vals = [16, 32, 64, 128];
%     L_vals = [9.25 18.5 37 74];
    T_vals = [.11, .14, .165, .17];
    T_dirs = {"T_.11", "T_.14", "T_.165", "T_.17"};
    dirs=sprintfc("/data/scc/thobi/211201_LongerTime/mxy_3.00/sqrtN_%d/",sqrtN_vals);
elseif dataset == "LinearTime"
    sqrtN_vals = [16, 32, 64, 128, 256];
%     L_vals = [9.25 18.5 37 74 148];
    T_vals = [.14, .17, .19, .21];
    T_dirs = {"T_.14", "T_.17", "T_.19", "T_.20"};
    dirs=sprintfc("/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_%d/",sqrtN_vals);
end
dirs_logtime=sprintfc("/data/scc/thobi/210727_LepriRuffo_GoodResolution/mxy_3.00/sqrtN_%d/",sqrtN_vals);    
% T_dirs=sprintfc("T_%.3g",T_vals); 

if t_scaling == "t"
    t_scale=@(t,L) t;
    x_str = "$t$";
elseif t_scaling == "tbyL"
    t_scale=@(t,L) t/L;
    x_str = "$t/L$";
elseif t_scaling == "tbyL2"
    t_scale=@(t,L) t/L^2;
    x_str = "$t/L^2$";
end



L_vals = [9.25, 18.5, 37, 74, 148];
fig_indices = [1 2 4 5];
sampfilename = "samp_Dynamics_M_TimeEvolution";
sampfilename_logtime = "samp_LepriRuffo_M_TimeEvolution";
N_T = numel(T_vals);
N_N = numel(sqrtN_vals);
figure
for i_T = 1:N_T
    T = T_vals(i_T);
    c_map=linspecer(N_N);
    absM_vec=zeros(1,N_N);
    for i_N = 1:N_N
        curdir = sprintf("%s/%s",dirs{i_N},T_dirs{i_T});
        curdir_logtime = sprintf("%s/%s",dirs_logtime{i_N},T_dirs{i_T});
        load(sprintf('%s/%s',curdir_logtime,sampfilename_logtime),...
            "absM_av");
        absM_vec(i_N) = absM_av;
%         absM_vec(i_N) = mean(mean(sqrt(Mx_collect.^2 + My_collect.^2)));
    end
    eta_fitob=fit_eta_Magnetization_FS(absM_vec(:),L_vals(1:N_N));
    eta = eta_fitob.eta;
    for i_N = 1:N_N
        sqrtN = sqrtN_vals(i_N);
        L = L_vals(i_N);
        curdir = sprintf("%s/%s",dirs{i_N},T_dirs{i_T});
        curdir_logtime = sprintf("%s/%s",dirs_logtime{i_N},T_dirs{i_T});
%         load(sprintf('%s/%s',curdir,sampfilename),...
%             "averaging_times","M_ang_MSD");%,"Mx_collect","My_collect");
%         t_lintime=averaging_times;
%         y_lintime=M_ang_MSD;
        load(sprintf('%s/%s',curdir_logtime,sampfilename_logtime),...
            "averaging_times","M_ang_MSD");
        t_logtime=averaging_times;
        y_logtime=M_ang_MSD;

        t = t_logtime;
        y = y_logtime;
%         [t,y]=combine_cf(t_logtime,y_logtime,t_lintime,y_lintime);
%         absM_av = mean(mean(sqrt(Mx_collect.^2 + My_collect.^2)));
        t_max = max(averaging_times);
        n_t = numel(averaging_times);
        t_max_exp = round(log10(t_max));
        t_max_noexp = t_max / 10^t_max_exp;
%         subfig_index = (i_N - 1) *N_T + i_T;
%         subplot(N_N,N_T,subfig_index)
        subfig_index = fig_indices(i_T);
        subplot(2,3,subfig_index)
        dispname=sprintf('$N = (%d)^2$',sqrtN);
        loglog(t,y*L^(2 - eta)./t.^2,...
            'Color',c_map(i_N,:),...
            'Displayname',dispname,...
            'LineWidth',2);
        hold on;

    
    
        
        set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'XGrid', 'off', 'YGrid', 'off', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
            'XTick',logspace(-10,4,15),...
            'LineWidth', .5,...
            'XScale','log','YScale','lin')

    end 
    h_axis = gca;
    h_title = title(sprintf("$T = %.3f$",T),'interpreter','latex');
    hXLabel = xlabel('t','interpreter','latex');
    hYLabel = ylabel('$L^{2-\eta} G(t)/t^2$','interpreter','latex');
    
    set(h_axis, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    
    set(h_title,'FontName','cmr12','fontsize',fontsize_labels)
    xlim([2e-2 1e2]);
    ylim([0,Inf])
%     xx=logspace(0,1);
%     plot(t_scale(xx,L),1e-6*xx.^2,'--',...
%         'LineWidth',2,...
%         'Color','black');
%     xx=logspace(2.5,4);
%     plot(t_scale(xx,L),1e-6*xx,'--',...
%         'LineWidth',2,...
%         'Color','black');
end

subplot (2, 3, [3 6]) % merge remaining subplots and put legend here
for i_N = 1:N_N
    dispname=sprintf('$N = (%d)^2$', sqrtN_vals(i_N));
    plot(nan, nan,'Color',c_map(i_N,:),'DisplayName',dispname,...
        'LineWidth',2) % plot nans (hack to generate correct legend but plot no data)
    hold on;
end
% hLegend = legend(sprintfc('$N = (%d)^2$', sqrtN_vals), 'Location', 'west',...
%     'interpreter','latex'); 
axis off
hLegend = legend('interpreter','latex','Location','West');
set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)

set(gcf,'units','centimeters','OuterPosition',[0 0 pagewidth_cm .5*pageheight_cm]);


figname=sprintf('%s/mxy_M_ang_MSD_yscaled_shorttime_%s',basedir,dataset);

fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

