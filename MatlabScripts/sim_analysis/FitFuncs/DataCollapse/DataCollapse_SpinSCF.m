clear all
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;

dataset_id='eq_mxy';
% dataset_id='dynamics_mxy';
% dataset_id='dynamics_mxy_fullT';
% dataset_id='dynamics_mxy_better_q';
% if (strcmp(dataset_id,'eq_mxy'))
%     mxydata=load('mxy/rho_3.00_eq.mat');
% elseif (strcmp(dataset_id,'dynamics_mxy'))
%     mxydata=load('mxy/rho_3.00_dynamics.mat');
% elseif (strcmp(dataset_id,'dynamics_mxy_fullT'))
%     mxydata=load('mxy/rho_3.00_dynamics_fullT.mat');
% elseif (strcmp(dataset_id,'dynamics_mxy_better_q'))
%     mxydata=load('mxy/rho_3.00_dynamics_better_q.mat');
% end
mxydata_name='mxy/rho_3.00_dynamics_LinearTime.mat';
xydata_name='xy/xy_dynamics_LinearTime.mat';
fmxydata_name='fmxy/fmxy_dynamics_LinearTime.mat';

mxyfit_name='mxy/rho_3.00_CritExpFit.mat';
xyfit_name='xy/xy_CritExpFit.mat';
fmxyfit_name='fmxy/fmxy_CritExpFit.mat';



for i_model = 3:3

    % figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=load(mxydata_name);
        fitdata=load(mxyfit_name);
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        rbin = data.('rbin');
        SCF_Spin_av = data.('SCF_Spin_av');
        storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_DataCollapse_SCF.mat';

        N_N = numel(sqrtN_vals);
        N_T = numel(T_vals);

        T_select=1:numel(T_vals);
        
        L_vals=[9.25,18.5,37,74,148];
        %% One should make sure to ignore FS effects
        for i_N = 1:N_N
            for i_T = 1:N_T
                index_max=find(rbin{i_N,i_T}>.35*L_vals(i_N),1);
                rbin{i_N,i_T} = rbin{i_N,i_T}(1:index_max);
                SCF_Spin_av{i_N,i_T} = SCF_Spin_av{i_N,i_T}(1:index_max);
            end
        end

    elseif (i_model == 2)
        curmodel="xy";
        curtitle="SXY model";
        
        data=load(xydata_name);
        fitdata=load(xyfit_name);
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        rbin = data.('rbin');
        SCF_Spin_av = data.('SCF_Spin_av');
        storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/xy/xy_DataCollapse_SCF.mat';

        N_N = numel(sqrtN_vals);
        N_T = numel(T_vals);

        T_select=1:numel(T_vals);
        
        L_vals=sqrtN_vals;
        % Bad last entries of SCF_Spin, trim that. Also, there are not
        % particles at all distances, so we ignore zero values in
        % SCF_Spin_av
        for i_N = 1:N_N
            for i_T = 1:N_T
                index_max=min(numel(rbin{i_N,i_T})-2, find(rbin{i_N,i_T}>.35*L_vals(i_N),1));
                rbin{i_N,i_T} = rbin{i_N,i_T}(1:index_max);
                SCF_Spin_av{i_N,i_T} = SCF_Spin_av{i_N,i_T}(1:index_max);
                
                nonzero_indices=find(SCF_Spin_av{i_N,i_T} ~= 0);
                rbin{i_N,i_T} = rbin{i_N,i_T}(nonzero_indices);
                SCF_Spin_av{i_N,i_T} = SCF_Spin_av{i_N,i_T}(nonzero_indices);
            end
        end
    elseif (i_model == 3)
        curmodel="fmxy";
        curtitle="FMXY model";
        
        data=load(fmxydata_name);
        fitdata=load(fmxyfit_name);
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        rbin = data.('rbin');
        SCF_Spin_av = data.('SCF_Spin_av');
%         if (strcmp(dataset_id,'eq_mxy'))
            storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/fmxy/fmxy_DataCollapse_SCF.mat';
%         end
        N_N = numel(sqrtN_vals);
        N_T = numel(T_vals);

        T_select=1:numel(T_vals);
        
        L_vals=[9.25,18.5,37,74,148];
        %% One should make sure to ignore FS effects
        for i_N = 1:N_N
            for i_T = 1:N_T
                index_max=find(rbin{i_N,i_T}>.35*L_vals(i_N),1);
                rbin{i_N,i_T} = rbin{i_N,i_T}(1:index_max);
                SCF_Spin_av{i_N,i_T} = SCF_Spin_av{i_N,i_T}(1:index_max);
            end
        end

    end
    eta_diff = .001;
    eta_max = 1;

%     [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);


    spline_curve=cell(N_N,N_T);
    spline_goodness=cell(N_N,N_T);
    spline_output=cell(N_N,N_T);
    spline_fiterr=zeros(N_N,N_T);
    
    
    for i_N = 1:length(sqrtN_vals)
        sqrtN_cur = sqrtN_vals(i_N);
        fprintf('++ Fitting. sqrtN = %d\n', sqrtN_cur);
        for i_T = 1:length(T_vals)
            T_cur = T_vals(i_T);
            r_cur = rbin{i_N,i_T};
            SCF_cur = SCF_Spin_av{i_N,i_T};
            [spline_curve{i_N,i_T}, spline_goodness{i_N,i_T}, spline_output{i_N,i_T}] = fit(r_cur(:),SCF_cur(:),'smoothingspline');
            spline_fiterr(i_N,i_T) = sqrt(sum((SCF_cur(:) - spline_curve{i_N,i_T}(r_cur)).^2/numel(r_cur)));
        end
    end
    
    
    eta_testvals = 0:eta_diff:eta_max;
    eta_vals = 0 * T_vals;
    eta_err = 0 * T_vals;
    error_vals = zeros(numel(T_vals),numel(eta_testvals));
    for i_T = 1:length(T_vals)
        T_cur = T_vals(i_T);
        fprintf('++ Collapsing. T = %.3f\n', T_cur);
%         error_vals = zeros(size(eta_testvals));
        curves=spline_curve(:,i_T);
        rbins=rbin(:,i_T);
        for i_eta = 1:numel(eta_testvals)
            eta=eta_testvals(i_eta);
            error_vals(i_T,i_eta) = CalculateCollapseError(curves,rbins,L_vals,1,eta,0);
        end
        [errmin,err_argmin] = min(error_vals(i_T,:));
        eta_vals(i_T) = eta_testvals(err_argmin);
        if (err_argmin ~= 1 && err_argmin ~= numel(eta_testvals))
            eta_err(i_T) = ((error_vals(i_T,err_argmin+1) + error_vals(i_T,err_argmin-1) - 2*error_vals(i_T,err_argmin) )...
                / (eta_diff^2))^(-1);
        end
%         figure(i_T)
%         plot(eta_vals,error_vals);
    end
    
    save(storefilename, 'spline_curve','spline_fiterr','spline_goodness',...
        'eta_vals','eta_err','error_vals','eta_testvals',...
        'sqrtN_vals','T_vals');
    
end