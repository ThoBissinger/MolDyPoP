%% Initialization
clear all
close all
initialization_script
basedir=sprintf('%s/plots/Spin_SCF',fig_base);
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;

% 
% mxydata=load('mxy/rho_3.00_eq.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat');
% xydata=load('xy/lf0_eq.mat'); xyfit=load('xy/lf0_CritExpFit.mat');

data_rootdir="/data/scc/thobi";
label='b';
lw=.8;
labelboxdim=[.05 .065];
boxbound=[.015 .02];
marksize=6; % marker size
for i_model = 1

    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        sqrtN_vals=[16, 32, 64, 128, 256];
        i_N_SCFfit=5;
        L_vals = 9.25/16*sqrtN_vals;
        symbols=['o','+','s','x','d'];
        data_dirs={"210715_LinearTimeSampling/mxy_3.00"};
        sampfilenames={"samp_Dynamics_collect"};
        simdir_helicity="210715_LinearTimeSampling/mxy_3.00";
        
        T_dirs = ["T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"];
        T_vals = [.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52 ];
        T_KT=.171;
        T_max_eta = .21;
        T_min_xi = .16;
        r_max_factor = .8;
        fit_confidentiality = .95;
        runmax=500;
        
        collapse_d_eta = .002;          % Data collapse sampling step
        collapse_interval_eta = .02;    % Data collapse sampling interval

%         T_dirs_eta = {"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23"};
%         T_vals_eta = [.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23];
%         T_dirs_xi = {"T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52" };
%         T_vals_xi = [.165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52 ];
    else
        curmodel="xy";
        curtitle="SXY model";
        
        r_min = 3;
        r_max = 55;
    end
%     N_T_eta = numel(T_vals_eta);
%     N_T_xi = numel(T_vals_xi);
    N_T = numel(T_vals);
    N_N = numel(sqrtN_vals);
    sampfile_helicity="samp_Dynamics_helicity";
    savefile_name=sprintf("data_%s_EtaFits",curmodel);

    %% 1 Gather data
    data_cf=cell(N_N,N_T);
    data_rbin=cell(N_N,N_T);
    data_absM=zeros(N_N,N_T);
    H_x_vals=zeros(N_N,N_T);
    H_y_vals=zeros(N_N,N_T);
    I_x_2_vals=zeros(N_N,N_T);
    I_y_2_vals=zeros(N_N,N_T);
    helicity_vals=zeros(N_N,N_T);
    for i_N = 1:N_N
        sqrtN=sqrtN_vals(i_N);
        fprintf('++ Loading data. sqrtN = %d\n', sqrtN);
        for i_T = 1:N_T
            T_dir=T_dirs(i_T);
            i_dir=1;
            curfile=sprintf('%s/%s/sqrtN_%d/%s/%s.mat',data_rootdir,data_dirs{i_dir},sqrtN,T_dir,sampfilenames{i_dir});
            while (~ isfile(curfile))
                i_dir = i_dir + 1;
                curfile=sprintf('%s/%s/sqrtN_%d/%s/%s.mat',data_rootdir,data_dirs{i_dir},sqrtN,T_dir,sampfilenames{i_dir});
            end
            load(curfile,'rbin_collect','SCF_Spin_av_collect','absM_av');
            data_rbin{i_N,i_T} = rbin_collect;
            data_cf{i_N,i_T} = SCF_Spin_av_collect;
            data_absM(i_N,i_T) = absM_av;
            
            filename_helicity=sprintf('%s/%s/sqrtN_%d/%s/samp_Dynamics_helicity.mat',...
                data_rootdir,simdir_helicity,sqrtN,T_dir);
            S_helicity=load(filename_helicity,...
                    "H_x","H_y","I_x_2","I_y_2");
            H_x_vals(i_N,i_T)=S_helicity.H_x;
            H_y_vals(i_N,i_T)=S_helicity.H_y;
            I_x_2_vals(i_N,i_T)=S_helicity.I_x_2; 
            I_y_2_vals(i_N,i_T)=S_helicity.I_y_2;
            helicity_vals(i_N,i_T)=1/2/L_vals(i_N)^2 * ...
                ( abs(H_x_vals(i_N,i_T) + H_y_vals(i_N,i_T)) - ...
                1/T_vals(i_T)*(I_x_2_vals(i_N,i_T) + I_y_2_vals(i_N,i_T)));
        end
    end
    
    %% Fitting
    T_vals_eta = T_vals(T_vals < T_max_eta);
    helicity_vals_eta=helicity_vals(:,T_vals <= max(T_vals_eta));
    N_T_eta = numel(T_vals_eta);
    eta_vals_FSMag=zeros(1,N_T_eta);
    eta_vals_FSMag_fitob=cell(1,N_T_eta);
    eta_vals_fitSCF=zeros(N_N,N_T_eta);
    eta_vals_fitSCF_fitob=cell(N_N,N_T_eta);
    eta_vals_fitSCF_ci=zeros(N_N,N_T_eta,2);
    eta_vals_fitSCF_sd=zeros(N_N,N_T_eta);

    eta_vals_FSSCF=zeros(1,N_T_eta);
    eta_vals_FSSCF_fitinterval=cell(N_N,N_T_eta);
    eta_vals_FSSCF_fiterr=cell(N_N,N_T_eta);

    spline_curve=cell(N_N,N_T_eta);
    spline_goodness=cell(N_N,N_T_eta);
    spline_output=cell(N_N,N_T_eta);
    spline_fiterr=zeros(N_N,N_T_eta);
    rbins=cell(1,N_N);
    curves=cell(1,N_N);
    
    % Eta fits
    for i_T = 1:N_T_eta
        T = T_vals_eta(i_T);
        fprintf('++ Fitting eta. T = %.3f\n', T);
        eta_vals_FSMag_fitob{i_T}=fit_eta_Magnetization_FS(data_absM(:,i_T),L_vals);
        eta_vals_FSMag(i_T)=eta_vals_FSMag_fitob{i_T}.eta;
        for i_N = 1:N_N
            r_vals = data_rbin{i_N,i_T};
            cf = data_cf{i_N,i_T};
            nonzero_ind = find(cf);
            r_vals = r_vals(nonzero_ind);
            cf = cf(nonzero_ind);
            
            fitob=fit_PowSCF(r_vals,cf,min(r_vals),r_max_factor * max(r_vals),false);
            eta_vals_fitSCF_fitob{i_N,i_T}=fitob;
            eta_vals_fitSCF(i_N,i_T)=fitob.eta;
            ci=confint(fitob,fit_confidentiality);
            eta_vals_fitSCF_ci(i_N,i_T,:)=ci(:,end);
            eta_vals_fitSCF_sd(i_N,i_T) = sqrt(runmax * (eta_vals_fitSCF_ci(i_N,i_T,2) - eta_vals_fitSCF_ci(i_N,i_T,1)).^2 / (2 * 1.96));

            [spline_curve{i_N,i_T}, spline_goodness{i_N,i_T}, spline_output{i_N,i_T}] = fit(r_vals(:),cf(:),'smoothingspline');
            spline_fiterr(i_N,i_T) = std(cf(:) - spline_curve{i_N,i_T}(r_vals)/sqrt(numel(r_vals)));

            rbins{i_N} = r_vals;
            curves{i_N} = spline_curve{i_N,i_T};
%             spline_fiterr(i_N,i_T) = sqrt(sum((cf(:) - spline_curve{i_N,i_T}(r_vals)).^2/numel(r_vals)));
        end
%         eta=eta_vals_FSMag(i_T);
        eta_best = round(100*eta_vals_FSMag(i_T))/100;
        eta_testvals = max(0,eta_best-collapse_interval_eta):collapse_d_eta:eta_best+collapse_interval_eta;
        eta_vals_FSSCF_fitinterval{i_T} = eta_testvals;
        error_vals = zeros(1,numel(eta_testvals));
    
        fprintf('++ Collapsing. T = %.3f\n', T);
%         error_vals = zeros(size(eta_testvals));
%         curves=spline_curve(:,i_T);
%         rbins=rbin(:,i_T);
        for i_eta = 1:numel(eta_testvals)
            eta=eta_testvals(i_eta);
            error_vals(i_T,i_eta) = CalculateCollapseError(curves,rbins,L_vals,1,eta,0);
        end
        [errmin,err_argmin] = min(error_vals(i_T,:));
        eta_vals_FSSCF(i_T) = eta_testvals(err_argmin);
%         if (err_argmin ~= 1 && err_argmin ~= numel(eta_testvals))
%             eta_err(i_T) = ((error_vals(i_T,err_argmin+1) + error_vals(i_T,err_argmin-1) - 2*error_vals(i_T,err_argmin) )...
%                 / (eta_diff^2))^(-1);
%         end
% %         figure(i_T)
% %         plot(eta_vals,error_vals);
    end
    xx=linspace(.16,.18,1e3);
    splinefit_FSMag = csapi(T_vals_eta(:),eta_vals_FSMag(:));
    yy=abs(fnval(splinefit_FSMag,xx)-.25);
    [minval,i_min] = min(abs(yy));
    eta_crit_FSMag = xx(i_min);

    splinefit_FSSCF = csapi(T_vals_eta(:),eta_vals_FSSCF(:));
    yy=abs(fnval(splinefit_FSSCF,xx)-.25);
    [minval,i_min] = min(abs(yy));
    eta_crit_FSSCF = xx(i_min);

    fprintf('Spline fit estimates of eta:\n');
    fprintf('  FS Mag: %.5f\n',eta_crit_FSMag);
    fprintf('  FS SCF: %.5f\n',eta_crit_FSSCF);
    fprintf('  FS SCF: ');
    for i_N = 1:N_N
        splinefit_SCF{i_N} = csapi(T_vals_eta(:),eta_vals_fitSCF(i_N_SCFfit,:));
        yy=abs(fnval(splinefit_SCF{i_N},xx)-.25);
        [minval,i_min] = min(abs(yy));
        eta_crit_SCF(i_N) = xx(i_min);
        fprintf('%.5f ',eta_crit_SCF(i_N));
    end
    fprintf('\n');

    %% Fitting of BKT-sqrt-fit
    fitfunc="1/4*(1-b/2/pi*sqrt(T_KT-x))";
    eta_f=@(b,T_c,x) .25*(1-b/2/pi*sqrt(T_c-x));
    ind=(T_vals_eta<1.05*T_KT) & (eta_vals_FSMag <= .25) .* (eta_vals_FSMag > .1);
    T_fit=T_vals_eta(ind);
    eta_fit=eta_vals_FSMag(ind);
    fitob_eta_BKT=fit(T_fit(:),eta_fit(:),fitfunc,...
        'StartPoint',[T_KT,1],"Lower",[max(T_fit),0],"Upper",[1.1*T_KT,Inf]);
    
    save(savefile_name);
    %%
    load(savefile_name);
    %% Plotting
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm],...
        'Resize','off');
    hold on
    c_map = linspecer(5);


    helicity_indices=find((abs(helicity_vals_eta(i_N_SCFfit,:)) < 1) .* ...
        (abs(helicity_vals_eta(i_N_SCFfit,:)) > 1e-5));
    y_h = T_vals_eta(helicity_indices)./(2*pi*helicity_vals_eta(i_N_SCFfit,helicity_indices));
    h_ind=find((abs(y_h < 1) .* ...
        abs(y_h) > 1e-5));
    heta_line_helicity = plot(T_vals_eta(h_ind),...
        y_h(h_ind),...
        'LineStyle', 'none', 'LineWidth', lw, ...
        'Color',c_map(1,:),...
        'Marker','d','MarkerSize',marksize,...
        'DisplayName','$1 / 2 \pi \beta \Upsilon$, Eqs. (7) \& (8)');


    heta_line_FS_Mag = plot(T_vals_eta,eta_vals_FSMag, ...
        'LineStyle', 'none', 'LineWidth', lw, ...
        'Color',c_map(2,:),...
        'Marker','^','MarkerSize',marksize,...
        'DisplayName','Fit to $\langle m\rangle$, eq. (22)');

    heta_line_SCF = plot(T_vals_eta,eta_vals_fitSCF(i_N_SCFfit,:),...
        'LineStyle', 'none', 'LineWidth', lw, ...
        'Color',c_map(3,:),...
        'Marker','o','MarkerSize',marksize,...
        'DisplayName','Fit to $C_m(r)$, Eq. (26)');
    
%     heta_line_SCF = errorbar(T_vals_eta,eta_vals_fitSCF(i_N_SCFfit,:),...
%         eta_vals_fitSCF_sd(i_N,:), ...
%         'LineStyle', 'none', 'LineWidth', lw, ...
%         'Color',c_map(3,:),...
%         'DisplayName','Fit to $C_m(r)$, Eq. (24)');
%     
    
%     heta_line_FS_SCF = plot(T_vals(find(eta_vals_FS_SCF)),eta_vals_FS_SCF(find(eta_vals_FS_SCF)), ...
    heta_line_FS_SCF = plot(T_vals_eta,eta_vals_FSSCF, ...
        'LineStyle', 'none', 'LineWidth', lw, ...
        'Color',c_map(4,:),...
        'Marker','v','MarkerSize',marksize,...
        'DisplayName','Collapse of $C_m(r,L)$, Eq. (27)');

    J_eff = T_vals_eta(2)./(2*pi*eta_vals_FSSCF(2));
    heta_line_lowT = plot(T_vals_eta,T_vals_eta./(2*pi*J_eff),...
        'LineStyle', ':', 'LineWidth', 1.01, ...
        'Color','k',...
        'DisplayName','$k_B T / 2 \pi J_{\textrm{eff}}$');

    TT=linspace(1.05*T_fit(1),fitob_eta_BKT.T_KT);
    heta_line_highTfit = plot(TT,eta_f(fitob.b,fitob.T_KT,TT),...
        'LineStyle', '-', 'LineWidth', 1.01, ...
        'Color','k',...
        'DisplayName','BKT fit');
    
    

    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', .5,'Xscale', 'lin','Yscale', 'lin')



    % Legend, axes etc
    hLegend = legend([heta_line_helicity,heta_line_FS_Mag,heta_line_SCF,heta_line_FS_SCF,heta_line_lowT],...
        'Location', 'northwest','interpreter','latex',...
        'NumColumns',1);
    xlim([0 T_max_eta]);
    ylim([0 .4]);

    hXLabel = xlabel('$T$','interpreter','latex');
    hYLabel = ylabel('$\eta$','interpreter','latex');

    % Font
    set(gca, 'FontName', 'cmr12','FontSize', fontsize_axis);
    set([hXLabel, hYLabel], 'FontName', 'cmr12','FontSize', fontsize_ax_labels);
    set(hLegend, 'FontName', 'cmr12','FontSize', fontsize_annotation)


    % Adjust axes properties
    
        
    
    if(saveswitch == 1)
        figname=sprintf('%s/%s_Eta_Compare',basedir,curmodel);
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end