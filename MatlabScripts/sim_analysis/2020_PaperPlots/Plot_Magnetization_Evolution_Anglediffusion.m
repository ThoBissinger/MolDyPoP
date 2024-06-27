clear;
close all;
initialization_script;

basedir=sprintf('%s/plots/Magnetization_Evolution',fig_base);
runmax=500;

dirs=["/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_16/",
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_32/",
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_64/",
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_128/",
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/"];
sqrtN_vals = [16, 32, 64, 128, 256];
T_dirs=["T_.07" "T_.11" "T_.14", "T_.17", "T_.18", "T_.19", "T_.21"];
T_vals = [.07, .11, .14, .17, .18, .19, .21];
% T_dirs=["T_.14", "T_.17", "T_.19", "T_.21"];
% T_vals = [.14, .17, .19, .21];

sampfilename = "sampling_output_Dynamics";
resultfilename = "samp_Dynamics_M_TimeEvolution";
N_T = numel(T_vals);
N_N = numel(sqrtN_vals);

MSD_max = pi^2/3;
fitfunc=@(c,t) c(1)*t .* (t <= 1e7*c(2)) + MSD_max * (t > 1e7*c(2));
D_coeff=zeros(N_N,N_T);
trans_time_coeff=zeros(N_N,N_T);
% finite_boundary_coeff=zeros(N_N,N_T);
coeffs=cell(N_N,N_T);
y_vals=cell(N_N,N_T);
t_vals=cell(N_N,N_T);
y_fit_vals=cell(N_N,N_T);
for i_T = 1:N_T
    for i_N = 1:N_N
        curdir = sprintf("%s/%s",dirs(i_N),T_dirs(i_T));
        load(sprintf('%s/%s',curdir,resultfilename),...
            "averaging_times","M_ang_MSD");
        t=averaging_times;
        y=M_ang_MSD;
        t_vals{i_N,i_T} = t;
        y_vals{i_N,i_T} = y;
        i_maxed = find(y > MSD_max,1);
        if isempty(i_maxed)
            D_initguess=y(end)/t(end);
        else
            D_initguess=y(ceil(3/4*i_maxed))/t(ceil(3/4*i_maxed));
        end
        trans_time_initguess = MSD_max / D_initguess;

        c0 = [y(20)/t(20), 1e-7*trans_time_initguess];
        [c,resnorm,~,exitflag,output] = lsqcurvefit(fitfunc,c0,t,y);
        D_coeff(i_N,i_T) = c(1);
        trans_time_coeff(i_N,i_T) = 1e7 * c(2);
        y_fit_vals{i_N,i_T} = fitfunc(c,t);
%         finite_boundary_coeff(i_N,i_T) = c(3);
        coeffs{i_N,i_T} = c;
            %"averaging_times","angles","Mx_collect","My_collect",...
            %"Mx_av","My_av","M_ang_MSD");
    end
end
return

basedir=sprintf('%s/plots/Magnetization_Evolution',fig_base);

dirs=["/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_16/",
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_32/",
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_64/",
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_128/",
    "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_256/"];
sqrtN_vals = [16, 32, 64, 128, 256];
T_dirs=["T_.14", "T_.17", "T_.19", "T_.21"];
T_vals = [.14, .17, .19, .21];
sampfilename = "sampling_output_Dynamics";
N_T = numel(T_vals);
N_N = numel(sqrtN_vals);
figure
for i_T = 1:N_T
    subplot(1,N_T,i_T)
    hold on;
    for i_N = 1:N_N
        curfile = sprintf("%s/%s/run_1/output/%s",dirs(i_N),T_dirs(i_T),sampfilename);
        sqrtN = sqrtN_vals(i_N);
        T = T_vals(i_T);
        load(curfile);
        t_max = max(averaging_times);
        t_max_exp = round(log10(t_max));
        t_max_noexp = t_max / 10^t_max_exp;
        plot(absM_av * cos(0:.01:2*pi),absM_av * sin(0:.01:2*pi),...
            'Color','red',...
            'LineWidth',2)
        
        plot(M(1,:),M(2,:),...
            'Color','blue',...
            'LineWidth',.3); 
        
        ax_max = 1.2;
        xlim([-ax_max,ax_max]);
        ylim([-ax_max,ax_max]);
        pbaspect([1 1 1])
    
    
        h_axis = gca;
        hXLabel = xlabel('$m_x$','interpreter','latex');
        hYLabel = ylabel('$m_y$','interpreter','latex');
    
        set([h_axis,hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_axis)
    
        set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'XGrid', 'off', 'YGrid', 'off', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
            'XTick', -1:.5:1, 'YTick', -1:.5:1, ...
            'LineWidth', .5)

        annotation_str = sprintf('$\\langle|m|\\rangle = %.2f$',absM_av);
        h_mval = text(.9, .1,annotation_str,'units','normalized',...
            'interpreter','latex',...
            'HorizontalAlignment','right',...
            'fontsize',fontsize_annotation,...
            'Color','red');
%         annotation_str = {sprintf('$\\langle|m|\\rangle = %.2f$',absM_av)};
%         dim=[.56 .13 .1 .2];
%         annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
%             'interpreter','latex', 'Units','normalized',...
%             'VerticalAlignment','middle', 'HorizontalAlignment','left',...
%             'Color','black','FontSize', fontsize_annotation,...
%             'BackgroundColor','white');
        if (i_N == 1)
            h_Tlabel = text(.5, 1.2,sprintf("$T = %.3f$",T)','units','normalized',...
                'interpreter','latex',...
                'HorizontalAlignment','center',...
                'fontsize',fontsize_labels);
        end
        if (i_T == 1)
            h_Nlabel = text(-.4, .5,sprintf("$N = (%d)^2$",sqrtN)','units','normalized',...
                'interpreter','latex',...
                'HorizontalAlignment','center',...
                'fontsize',fontsize_labels,...
                'Rotation',90);
%             set(h_Nlabel,'Rotation',90);
        end
        
    end 
end
% gcf

% annotation_str = {sprintf('$T = %.3f$',T),...
%     sprintf('$\\langle|m|\\rangle = %.2f$',absM_av),...
%     sprintf('$t_{\\textrm{max}} = %.1f \\cdot 10^{%d}$',t_max_noexp,t_max_exp)};
% dim=[.56 .13 .1 .2];
% annotation('textbox',dim, 'String', annotation_str, 'FitBoxToText','on',...
%     'interpreter','latex',...
%     'VerticalAlignment','middle', 'HorizontalAlignment','left',...
%     'Color','black','FontSize', fontsize_annotation,...
%     'BackgroundColor','white');
set(gcf,'units','centimeters','OuterPosition',[0 0 pagewidth_cm pageheight_cm]);


%     figname=sprintf('%s/mxy_sqrtN_%d_T_%.3f',basedir,sqrtN,T);
figname=sprintf('%s/mxy_M_Evolution',basedir);
fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

