clear
curdir="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/SCF_FS";
addpath("/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/FitFuncs/DataCollapse");
cd(curdir);
run ../initialization_script
saveswitch=1;
basedir=sprintf('%s/Plots',curdir);

model = "mxy";
if model == "mxy"
    sqrtN_vals=[16, 32, 64, 128];
    L_vals=9.25*2.^[0:3];
    T_str=[".11" ".13" ".14" ".15" ".16" ".17" ".21"];
    T_vals=str2double(T_str);
    filedir="/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00";
    filename="samp_Dynamics_grSCF.mat";
end
i_N=4;
sqrtN=sqrtN_vals(i_N);
L=L_vals(i_N);
N_T = numel(T_vals);
N_N = numel(sqrtN_vals);
pathbase=sprintf("%s/Plots",curdir);
%% Load data
Cm_cell=cell(N_N,N_T);
rbin_cell=cell(N_N,N_T);
gr_cell=cell(N_N,N_T);
fitob_cell=cell(N_N,N_T);
absM_vals=zeros(N_N,N_T);
eta_vals=zeros(N_N,N_T);
eta_M_vals=zeros(1,N_T);
eta_FS_vals=zeros(1,N_T);
for i_N = 1:N_N
    for i_T = 1:N_T
        curfile=sprintf("%s/sqrtN_%d/T_%s/%s",filedir,sqrtN_vals(i_N),T_str(i_T),filename);
        S=load(curfile);
        Cm_cell{i_N,i_T}=S.Cm;
        rbin_cell{i_N,i_T}=S.rbin;
        gr_cell{i_N,i_T}=S.gr;
        absM_vals(i_N,i_T)=S.absM;
        fitob=fit(log(S.rbin(S.gr>.4))',log(S.Cm(S.gr>.4)./S.gr(S.gr>.4))',"a-eta*x", ...
            "StartPoint",[0, .25],"Lower",[-Inf,0]);
        eta_vals(i_N,i_T)=fitob.eta;
        fitob_cell{i_N,i_T}=fitob;
    end
end
for i_T = 1:N_T
    fitob_M=fit(log(2*sqrtN_vals.^2)',log(absM_vals(:,i_T)),"a-eta*x/4", ...
            "StartPoint",[0, .25],"Lower",[-Inf,0]);
    eta_M_vals(i_T)=fitob_M.eta;
end

%% Data Collapse
eta_diff = .002;
eta_max = .5;



spline_curve=cell(N_N,N_T);
spline_goodness=cell(N_N,N_T);
spline_output=cell(N_N,N_T);
spline_fiterr=zeros(N_N,N_T);


for i_N = 1:N_N
    sqrtN_cur = sqrtN_vals(i_N);
    fprintf('++ Fitting. sqrtN = %d\n', sqrtN_cur);
    for i_T = 1:N_T
        T_cur = T_vals(i_T);
        r_cur = rbin_cell{i_N,i_T};
        SCF_cur = Cm_cell{i_N,i_T};
        [spline_curve{i_N,i_T}, spline_goodness{i_N,i_T}, spline_output{i_N,i_T}] = fit(r_cur(:),SCF_cur(:),'smoothingspline');
        spline_fiterr(i_N,i_T) = sqrt(sum((SCF_cur(:) - spline_curve{i_N,i_T}(r_cur)).^2/numel(r_cur)));
    end
end


% eta_testvals = 0:eta_diff:eta_max;
eta_FS_vals = 0 * T_vals;
% eta_err = 0 * T_vals;
% error_vals = zeros(numel(T_vals),numel(eta_testvals));
for i_T = 1:length(T_vals)
    eta_cur=eta_M_vals(i_T);
%     eta_testvals=eta_cur/2:eta_diff:2*eta_cur;
%     error_vals = zeros(numel(T_vals),numel(eta_testvals));
    T_cur = T_vals(i_T);
    fprintf('++ Collapsing. T = %.3f\n', T_cur);
%         error_vals = zeros(size(eta_testvals));
    curves=spline_curve(:,i_T);
    rbins=rbin_cell(:,i_T);
    eta_FS_vals(i_T)=fminbnd(@(eta) CalculateCollapseError(curves,rbins,L_vals,1,eta,0),0,1);
%     for i_eta = 1:numel(eta_testvals)
%         eta=eta_testvals(i_eta);
%         error_vals(i_T,i_eta) = CalculateCollapseError(curves,rbins,L_vals,1,eta,0);
%     end
%     [errmin,err_argmin] = min(error_vals(i_T,:));
%     eta_vals(i_T) = eta_testvals(err_argmin);
%     if (err_argmin ~= 1 && err_argmin ~= numel(eta_testvals))
%         eta_err(i_T) = ((error_vals(i_T,err_argmin+1) + error_vals(i_T,err_argmin-1) - 2*error_vals(i_T,err_argmin) )...
%             / (eta_diff^2))^(-1);
%     end
%         figure(i_T)
%         plot(eta_vals,error_vals);
end
%% Create basic plot
figure
hold on
N_sqrtN = numel(sqrtN_vals);
c_map = linspecer(numel(sqrtN_vals));
for i_T = 1:N_T
%     i_T = T_select(i);
    T = T_vals(i_T);
    eta=eta_FS_vals(i_T);
    figure
    for i_N = 1:numel(sqrtN_vals)
        sqrtN = sqrtN_vals(i_N);
        SCF=Cm_cell{i_N,i_T}./gr_cell{i_N,i_T};
        
        r_vals = rbin_cell{i_N,i_T};
        r_min = r_vals(2);
        r_max = r_vals(end);

        rr = linspace(r_vals(1),r_vals(end));
%             curpowfit=fit_PowSCF(r_vals,SCF,r_min,r_max,0);
%             powfit{i} = curpowfit;
%                 hFit_line(i_N) = line(rr,curpowfit(rr));
%                 set(hFit_line(i_N), 'LineStyle', '--')
%             set(hFit_line(i_N), ...
%                 'LineWidth', 2, ...
%                 'HandleVisibility','Off',...
%                 'Color',c_map(i_N,:));

        dispname=sprintf('$N = %d^2$', sqrtN);
        hSCF_line(i_N) = line(r_vals(1:end)/L_vals(i_N),SCF(1:end)*L_vals(i_N)^eta);
        set(hSCF_line(i_N), ...
            'LineStyle', 'none', 'LineWidth', 1, ...
            'DisplayName',dispname,...
            'Marker', '.', 'MarkerSize', 10, ...
            'Color',c_map(i_N,:),'MarkerEdgeColor', c_map(i_N,:),'MarkerFaceColor', c_map(i_N,:))
    end






    % Legend, axes etc
    hLegend = legend('Location', 'Southwest','interpreter','latex','FontSize', 20,...
        'NumColumns',1);
    xlim([1e-2 Inf])
%         ylim([1e-2 3e0]);
    % hlegend(2).LineStyle = '-';
    % lineEntry = findobj(hLegend.EntryContainer, 'Object',hMag_dots(2));
    % entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    % entryMarker.LineWidth = .5;
    
    hXLabel = xlabel('$r/L$','interpreter','latex');
    hYLabel = ylabel('$L^\eta C_m(r)$','interpreter','latex');
%         hTitle = title(curtitle);

    % Font
    set(gca, 'FontName', 'cmr12')
    set([hXLabel, hYLabel], 'FontName', 'cmr12')
    set(gca, 'FontSize', 8)
    set(hLegend, 'FontSize', 12)
    set([hXLabel, hYLabel], 'FontSize', 12)
%         set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')


    % Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'log','Yscale', 'lin')
    ax_full = gca;
    
    text(.05,.5,sprintf('$T = %.3f$',T),...
        'Units','normalized',...
        'VerticalAlignment','top', 'HorizontalAlignment','left',...
        'interpreter','latex',...
        'Color','Red','FontSize', 10);
    text(.05,.45,sprintf('$\\eta = %.2f$',eta),...
        'Units','normalized',...
        'VerticalAlignment','top', 'HorizontalAlignment','left',...
        'interpreter','latex',...
        'Color','Red','FontSize', 10);

    figname=sprintf('%s/%s_SpinSCF_FS_T_%.3f',pathbase,model,T);
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
end
