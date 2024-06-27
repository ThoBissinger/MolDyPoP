pathbase="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/Static/Magnetization";
cd(pathbase)
run('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/initialization_script.m')

model="xy";

savefile=sprintf("data_%s",model);
if model == "mxy"
    sqrtN_vals = [16,32,64,128,256];
    T_str = [".01" ".03" ".05" ".07" ".09" ".11" ".13" ".14" ".15" ".155" ".16" ".165" ".17" ".175" ".18" ".185" ".19" ".195" ".20" ".205" ".21" ".22" ".23" ".24" ".25" ".27" ".29" ".31" ".33" ".35" ".37" ".40" ".43" ".46" ".49" ".52" ];
    filebase = "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_";
    fileend = "/samp_Dynamics";
    T_KT=.173;
elseif model == "xy"
    sqrtN_vals = [16,32,64,128];
%     T_str = [".10" ".20" ".30" ".40" ".50" ".60" ".70" ".80" ".85" ".91" ".95" "1.00"];
    T_str = [".10" ".20" ".30" ".40" ".50" ".60" ".70" ".80" ".85" ".87" ".89" ".90" ".91" ".93" ".95" ".97" "1.00" "1.03" "1.06" "1.09" "1.10" "1.12" "1.15" "1.18" "1.20" "1.21" "1.24" "1.30" "1.40" "1.50" "1.60" "1.70" "1.80" "1.90" "2.00"];
%     filebase = "/data/scc/thobi/220201_ReducedSmapleStepDeltat/xy_s/sqrtN_";
    filebase = "/data/scc/thobi/201207_equilibration/xy_s/scale/sqrtN_";
    filebase_anneal = "/data/scc/thobi/201207_equilibration/xy_s/anneal/sqrtN_";
    fileend = "/samp_eq";
    T_KT=.89;
end
linewidth=1.3;
T_vals=str2double(T_str);
N_N=numel(sqrtN_vals);
N_T=numel(T_vals);

% filebase = "/data/scc/thobi/210715_LinearTimeSampling//mxy_3.00/sqrtN_";
filemid = "/T_";

% absM_vals=cell2mat(data.('absM_av'));
% absM_var_vals=cell2mat(data.('absM_var'));
% eta_vals_FS_Mag = FS_Mag_fitdata.('eta_vals');

absM_vals=zeros(N_N,N_T);
absM_var_vals=zeros(N_N,N_T);
eta_vals_Maglog = zeros(N_N,N_T);
eta_vals_FS_Mag = zeros(1,N_T);
a_vals_FS_Mag = zeros(1,N_T);

for i_N = 1:N_N
    
    sqrtN=sqrtN_vals(i_N);
    fprintf("%d ",sqrtN)
    for i_T = 1:N_T
        T = T_vals(i_T);
        fprintf("%.3f ",T)
        if model == "xy" && sqrtN == 16
            filenames(i_N,i_T) = filebase_anneal + sqrtN + filemid + T_str(i_T) + fileend;
        else
            filenames(i_N,i_T) = filebase + sqrtN + filemid + T_str(i_T) + fileend;
        end
        S{i_N,i_T}=load(filenames(i_N,i_T),"absM_av","absM_var");

        absM_vals(i_N,i_T)=S{i_N,i_T}.absM_av;
        absM_var_vals(i_N,i_T)=S{i_N,i_T}.absM_var;
        eta_vals_Maglog(i_N,i_T)=-4*log(absM_vals(i_N,i_T))/log(2*sqrtN^2);
    end
    fprintf("\n");
end
for i_T = 1:N_T
    fitob=fit(sqrtN_vals(:),absM_vals(:,i_T),"a*x^(-eta/2)","StartPoint",[1,.25]);
    eta_vals_FS_Mag(i_T)=fitob.eta;
    a_vals_FS_Mag(i_T)=fitob.a;
end


save(savefile,"absM_vals","absM_var_vals","eta_vals_Maglog",...
    "eta_vals_FS_Mag","a_vals_FS_Mag");
return
for i_model = [1]

    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        runmax=500;
        
        data=matfile(mxydata_name);
        FS_Mag_fitdata=matfile(mxyfit_FSMag_name);
        
        L_vals=[9.25,18.5,37,74,148];
        T_max = .4;
        FSplot_min = .15;
        FSplot_max = .255;
        FS_Tstep = .025;
        T_offset = .0025; % For text in inset
        inset_Tmin = .138;
        inset_Tmax = .252;
        
        inset_T_ticks=[.15 .19 .23];
    elseif (i_model == 2)
        curmodel="xy";
        curtitle="SXY model";
        runmax=250;
        
        data=matfile(xydata_name);
        FS_Mag_fitdata=matfile(xyfit_FSMag_name);
%         T_vals=data.('T_vals');
%         sqrtN_vals=data.('sqrtN_vals');
%         absM_av=cell2mat(data.('absM_av'));
%         absM_var=cell2mat(data.('absM_var'));
        L_vals=[16 32 64 128 256];
        T_max = 3.5;
        FSplot_min = 1.2;
        FSplot_max = 2.1;
        FS_Tstep = .5;
        T_offset = .03; % For text in inset
        inset_Tmin = 1.3;
        inset_Tmax = 2.1;

        inset_T_ticks=[1.4 1.7 2.0];
    elseif (i_model == 3)
        curmodel="fmxy";
        curtitle="FMXY model";
        runmax=500;
        
        data=matfile(fmxydata_name);
        FS_Mag_fitdata=matfile(fmxyfit_FSMag_name);
        
%         T_vals=data.('T_vals');
%         sqrtN_vals=data.('sqrtN_vals');
%         absM_av=cell2mat(data.('absM_av'));
%         absM_var=cell2mat(data.('absM_var'));
        L_vals=[9.25,18.5,37,74,148];
        T_max = .4;
        FSplot_min = .15;
        FSplot_max = .255;
        FS_Tstep = .025;
        T_offset = .0025; % For text in inset
        inset_Tmin = .138;
        inset_Tmax = .252;
        inset_T_ticks=[.15 .19 .23];
    elseif (i_model == 4)
        curmodel="xy_s";
        curtitle="SXY model";
        runmax=250;
        
        data=matfile(xysdata_name);
        FS_Mag_fitdata=matfile(xysfit_FSMag_name);
%         T_vals=data.('T_vals');
%         sqrtN_vals=data.('sqrtN_vals');
%         absM_av=cell2mat(data.('absM_av'));
%         absM_var=cell2mat(data.('absM_var'));
        L_vals=[16 32 64 128 256];
        T_max = 2;
        FSplot_min = 1.2;
        FSplot_max = 2.1;
        FS_Tstep = .1;
        T_offset = .03; % For text in inset
        inset_Tmin = .78;
        inset_Tmax = 1.3;
        inset_T_ticks=[.8 1 1.2];
        
    end
    inset_N_ticks=[1e2 1e3 1e4];
    absM_vals=cell2mat(data.('absM_av'));
    absM_var_vals=cell2mat(data.('absM_var'));
    eta_vals_FS_Mag = FS_Mag_fitdata.('eta_vals');
    labels=["$N = (16)^{2}$", "$N = (32)^{2}$", "$N = (64)^{2}$", "$N = (128)^{2}$", "$N = (256)^{2}$"];

    [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_vals,sqrtN_vals,sqrt(2));
    %% 1 Create basic plot
    figure
    hold on
    N_sqrtN = numel(sqrtN_vals);
    c_map = linspecer(N_sqrtN);
    for i_N = 1 : N_sqrtN
        fitfunc=mag_fitfuncs{i_N};

        TT=linspace(.7*T_star(i_N), T_C(i_N));
        hfit_line(i_N) = line(TT,fitfunc(T_C(i_N),TT));
        set(hfit_line(i_N), ...
            'LineStyle', '-', 'LineWidth', 1.2, ...
            'HandleVisibility','off',...
            'Color', c_map(i_N,:))

        
    end
    for i_N = 1 : N_sqrtN
        hMag_dots(i_N) = line(T_vals,absM_vals(i_N,:));
%         set(hMag_dots(i_N), ...
%             'LineStyle', 'none', ...
%             'LineWidth', .4, ...
%             'Marker', 'o', 'MarkerSize', 4, ...
%             'DisplayName', labels(i_N), ...
%             'Color',c_map(i_N,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))
        set(hMag_dots(i_N), ...
            'LineStyle', 'none', ...
            'LineWidth', .8,...
            'Marker', 's', 'MarkerSize', 4, ...
            'DisplayName', labels(i_N), ...
            'Color',c_map(i_N,:))
        
    end
    for i_N = 1 : N_sqrtN
        hcrossover_line(i_N) = line(T_star(i_N),crossover_M(i_N));
        set(hcrossover_line(i_N), ...
            'LineStyle', 'none', ...
            'LineWidth', 1, ...
            'Marker', 'o', 'MarkerSize', 5.5, ...
            'HandleVisibility','off',...
            'Color', c_map(i_N,:), 'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))
    end

    hLegend = legend('Location', 'SouthWest','interpreter','latex');
    % Legend, axes etc
%         hLegend = legend([hMag_dots(1), hMag_dots(2), hMag_dots(3), hMag_dots(4), hMag_dots(5)], ...
%             '$N = (16)^2$', '$N = (32)^2$', '$N = (64)^2$', '$N = (128)^2$', '$N = (256)^2$', ...
%             'Location', 'SouthWest','interpreter','latex','FontSize', 12);
    % hlegend(2).LineStyle = '-';
    % lineEntry = findobj(hLegend.EntryContainer, 'Object',hMag_dots(2));
    % entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    % entryMarker.LineWidth = .5;
    hXLabel = xlabel('$T$','interpreter','latex');
    hYLabel = ylabel('$\langle m\rangle$','interpreter','latex');
%     hTitle = title(curtitle);
    xlim([0, T_max]);

    % Font
    set_fonts_default
%     set([, 'FontSize', 7)
%     set(hTitle, 'FontSize', 7, 'FontWeight' , 'bold')


    % Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:.2:1, ...
        'LineWidth', .5)
    ax_full = gca;
    
    % Inset
    ax_inset = axes('Position',[.6 .6 .28 .28]);
    xlabelstr=sprintfc("(%d)^2",sqrtN_vals);
    set(ax_inset, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'off', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'YTick', inset_T_ticks, 'XTick', 16^2*4.^[0:4],...
        'XTickLabel',xlabelstr,...
        'LineWidth', .5, 'Xscale', 'log')
    ylim([FSplot_min, FSplot_max]);
    
    hTC_line = line(sqrtN_vals.^2,T_C);
    hTstar_line = line(sqrtN_vals.^2,T_star);
    hTKT_line = line(sqrtN_vals.^2,T_KT);
    
    
    c_map = linspecer(3);
    set(hTC_line, ...
        'LineStyle', '--', 'LineWidth', .8,...
        'Marker', 's', 'MarkerSize', 4, ...
        'Color',c_map(1,:))
    set(hTstar_line, ...
        'LineStyle', '--', 'LineWidth', .8,...
        'Marker', 's', 'MarkerSize', 4, ...
        'Color',c_map(2,:))
    set(hTKT_line, ...
        'LineStyle', '--', 'LineWidth', .8,...
        'Marker', 's', 'MarkerSize', 4, ...
        'Color',c_map(3,:))
%     c_map = linspecer(3);
%     set(hTC_line, ...
%         'LineStyle', '--', 'LineWidth', 1.2, ...
%         'Marker', '^', 'MarkerSize', 5, ...
%         'Color',c_map(1,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(1,:))
%     set(hTstar_line, ...
%         'LineStyle', '--', 'LineWidth', 1.2, ...
%         'Marker', 's', 'MarkerSize', 5, ...
%         'Color',c_map(2,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(2,:))
%     set(hTKT_line, ...
%         'LineStyle', '--', 'LineWidth', 1.2, ...
%         'Marker', 'diamond', 'MarkerSize', 5, ...
%         'Color',c_map(3,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(3,:))
    
%     hXLabel = xlabel('$N$','interpreter','latex');
%     hYLabel = ylabel('$T$','interpreter','latex');
    set(gca, 'FontSize', fontsize_axis, 'FontName', 'cmr12');
%     set(hLegend, 'FontSize', fontsize_annotation, 'FontName', 'cmr12');
%     set_fonts_default
    
    % Text
    text(sqrtN_vals(2)^2,T_C(2)+T_offset,'$T_C$',...
        'VerticalAlignment','bottom','interpreter','latex',...
        'Color',c_map(1,:),'FontSize', fontsize_annotation);
    text(sqrtN_vals(2)^2,T_star(2)+T_offset,'$T^*$',...
        'VerticalAlignment','bottom','interpreter','latex',...
        'Color',c_map(2,:),'FontSize', fontsize_annotation);
    text(sqrtN_vals(2)^2,T_KT(2)-T_offset,'$T_{BKT}$',...
        'VerticalAlignment','top','interpreter','latex',...
        'Color',c_map(3,:),'FontSize', fontsize_annotation);
    
    
    ylim([inset_Tmin,inset_Tmax]);
    xlim([.9*sqrtN_vals(1)^2, 1.1*sqrtN_vals(end)^2]);
    
    NW = [min(xlim) max(ylim)];
    SE = [max(xlim) min(ylim)];
    SW = [min(xlim) min(ylim)];
    text(SE(1),SE(2),'$N$',...
        'VerticalAlignment','bottom','HorizontalAlignment','right',...
        'interpreter','latex',...
        'Color','black','FontSize', fontsize_axis);
    text(NW(1),NW(2),'\ $T$',...
        'VerticalAlignment','top','HorizontalAlignment','left',...
        'interpreter','latex',...
        'Color','black','FontSize', fontsize_axis);

%     pos = get(gca, 'position');
%     dim = [pos(3)-.05 .8*pos(2) pos(3) pos(4)]; %top right
% %     dim = [pos(1)+.025, .1, pos(3), .3*pos(4)]; % bottom left
% %     dim = [pos(3)-.05, .1, pos(3), .3*pos(4)]; % bottom right
% %     str = {dt_str,tmax_str,om_str,sig_str,tau_str};
%     str = {'\color{blue}$T_C$','\color{red}$T^*$','\color{green}$T_{\textrm{BKT}}$'};
%     annotation('textbox',dim,'String',str,'FitBoxToText','on',...
%         'interpreter','latex');

    
%     set(gcf,'units','points','OuterPosition',[0 0 240 3/4*240]);
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

    figname=sprintf('%s/%s_Magnetization',basedir,curmodel);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        print(sprintf('%s.eps',figname),'-depsc');
    end

    %% 2 Errors of magnetization
    figure
    hold on
    N_sqrtN = numel(sqrtN_vals);
    c_map = linspecer(N_sqrtN);
    for i_N = 1 : N_sqrtN
        hMag_dots(i_N) = line(T_vals,absM_var_vals(i_N,:)/sqrt(runmax));
        dispname=sprintf('$N=(%d)^2$',sqrtN_vals(i_N));
        set(hMag_dots(i_N), ...
            'LineStyle', '-', ...
            'Marker', 'o', 'MarkerSize', 6, ...
            'DisplayName', dispname, ...
            'LineWidth',1.5,...
            'Color',c_map(i_N,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))
    end


    % Legend, axes etc
        hLegend = legend('Location', 'NorthWest','interpreter','latex','FontSize', 20);
    % hlegend(2).LineStyle = '-';
    % lineEntry = findobj(hLegend.EntryContainer, 'Object',hMag_dots(2));
    % entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    % entryMarker.LineWidth = .5;
    hXLabel = xlabel('$T$','interpreter','latex');
    hYLabel = ylabel('$m_{\textrm{err}}$','interpreter','latex');
%     hTitle = title(curtitle);
    xlim([0, T_max]);

    % Font
    set_fonts_default


    % Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:1e-4:2e-3, ...
        'LineWidth', .5)
    ax_full = gca;

    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    
    figname=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/static/%s_Magnetization_errors',curmodel);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        print(sprintf('%s.eps',figname),'-depsc');
    end



    %% 3 Magnetization scaled with eta
    c_map = linspecer(N_sqrtN);
    figure
    
    for i_N = 1 : N_sqrtN
        
        hMag_dots(i_N) = line(T_vals,absM_vals(i_N,:).*((2*L_vals(i_N)).^(.5*eta_vals_FS_Mag)));
        set(hMag_dots(i_N), ...
            'LineStyle', '--', ...
            'Marker', 'o', 'MarkerSize', 7, ...
            'LineWidth', 2, ...
            'DisplayName', labels(i_N), ...
            'Color',c_map(i_N,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:)) 
    end
    hLegend = legend('Location', 'NorthWest','interpreter','latex','FontSize', 12);

    hXLabel = xlabel('$T$','interpreter','latex');
    hYLabel = ylabel('$m (2L)^{\eta/2}$','interpreter','latex');
%     hTitle = title(curtitle);
    xlim([0, T_max]);

    % Font
    set_fonts_default


    % Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'of', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'LineWidth', .5)
    ax_full = gca;

    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    
    figname=sprintf('%s/%s_Magnetization_Collapse',basedir,curmodel);
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        print(sprintf('%s.eps',figname),'-depsc');
    end
end
    



return
%% 4 Comparing Magnetizations for frozen and mobile
figure(1e4)
load(mxydata_name,'T_vals');
load(mxydata_name,'sqrtN_vals');
load(mxydata_name,'absM_vals');
absM_av_mxy=cell2mat(absM_vals);
load(fmxydata_name,'absM_vals');
absM_av_fmxy=cell2mat(absM_vals);
L_vals=[9.25,18.5,37,74,148];
T_max = .4;
labels=['N = (16)^{2}', 'N = (32)^{2}', 'N = (64)^{2}', 'N = (128)^{2}', 'N = (256)^{2}'];
hold on
N_sqrtN = numel(sqrtN_vals);
c_map = linspecer(N_sqrtN);
for i_N = 1 : N_sqrtN
    hMag_dots_mxy(i_N) = line(T_vals,absM_av_mxy(i_N,:) - absM_av_fmxy(i_N,:));
    dispname=sprintf('$N=(%d)^2$',sqrtN_vals(i_N));
    set(hMag_dots_mxy(i_N), ...
        'LineStyle', '--', ...
        'Marker', 'o', 'MarkerSize', 4, ...
        'DisplayName', dispname, ...
        'LineWidth',2,...
        'Color',c_map(i_N,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))
end


% Legend, axes etc
    hLegend = legend('NumColumns',1,...
        'Location', 'NorthWest','interpreter','latex','FontSize', 12);
% hlegend(2).LineStyle = '-';
% lineEntry = findobj(hLegend.EntryContainer, 'Object',hMag_dots(2));
% entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
% entryMarker.LineWidth = .5;
hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$\langle m\rangle$','interpreter','latex');
xlim([0, T_max]);

% Font
set_fonts_default

% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', -1e-2:1e-3:1e-2, ...
    'LineWidth', .5)
ax_full = gca;

set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

% Saving data
% figname='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/static/fmxy_mxy_MagnetizationCompare';
figname=sprintf('%s/fmxy_mxy_MagnetizationCompare',basedir);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

%% Figure comparing different exponents
fitfunc{1}="-a/2*log(b*x)";
fitfunc{2}="-a/2*log(1.3585*x)";
fitfunc{3}="-a/2*log(sqrt(2)*x)";
fitfunc{4}="-a/2*log(2*x)";
fitfunc{5}="-a/2*log(x/2)";
StartPoints={[.25,1.3585],.25,.25,.25,.25};
fitobs=cell(4,numel(T_vals));
a_vals=zeros(4,numel(T_vals));
b_vals=zeros(1,numel(T_vals));
for i_T = 1:numel(T_vals)
    fprintf('%d ',i_T);
    for i = 1:numel(fitfunc)
        fitobs{i,i_T} = fit(sqrtN_vals(:),log(absM_vals(:,i_T)),fitfunc{i},'StartPoint',StartPoints{i});
        a_vals(i,i_T) = fitobs{i,i_T}.a;
    end
    b_vals(i_T) = fitobs{1,i_T}.b
end
