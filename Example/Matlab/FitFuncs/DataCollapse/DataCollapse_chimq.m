clear all
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;

% dataset_id='eq_mxy';
dataset_id='dynamics_mxy';
dataset_id='dynamics_LinearTime_mxy';
% dataset_id='dynamics_mxy_fullT';
% dataset_id='dynamics_mxy_better_q';
if (strcmp(dataset_id,'eq_mxy'))
    mxydata=load('mxy/rho_3.00_eq.mat');
elseif (strcmp(dataset_id,'dynamics_mxy'))
    mxydata=load('mxy/rho_3.00_dynamics.mat');
elseif (strcmp(dataset_id,'dynamics_LinearTime_mxy'))
    mxydata=load('mxy/rho_3.00_dynamics_LinearTime.mat');
elseif (strcmp(dataset_id,'dynamics_mxy_fullT'))
    mxydata=load('mxy/rho_3.00_dynamics_fullT.mat');
elseif (strcmp(dataset_id,'dynamics_mxy_better_q'))
    mxydata=load('mxy/rho_3.00_dynamics_better_q.mat');
end
% mxydata=load('mxy/rho_3.00_integ.mat');
% mxydata=load('mxy/rho_3.00_dynamics.mat');
mxyfit=load('mxy/rho_3.00_CritExpFit.mat');

xydata=load('xy/xy_dynamics_LinearTime.mat');
% xydata=load('xy/lf0_qreduced.mat');
xyfit=load('xy/xy_CritExpFit.mat');


for i_model = 2:2

    % figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        data=mxydata;
        fitdata=mxyfit;
        
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        qbin = data.('qbin');
        chimxq_av = data.('chimxq_av');
        chimyq_av = data.('chimxq_av');
        chimparq_av = data.('chimxq_av');
        chimperpq_av = data.('chimxq_av');
%         if (strcmp(dataset_id,'eq_mxy'))
            storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_DataCollapse_chimq.mat';
%         end
        N_N = numel(sqrtN_vals);
        N_T = numel(T_vals);

        T_select=1:numel(T_vals);
        
        L_vals=[9.25,18.5,37,74,148];
%         L_vals=[9.25,18.5,37,74,148];


    else
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % TODO
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        curmodel="xy";
        curtitle="SXY model";
        
        data=xydata;
        fitdata=xyfit;
        
        T_vals=data.('T_vals');
        sqrtN_vals=data.('sqrtN_vals');
        qbin = data.('qbin');
        chimxq_av = data.('chimxq_av');
        chimyq_av = data.('chimxq_av');
        chimparq_av = data.('chimxq_av');
        chimperpq_av = data.('chimxq_av');
%         if (strcmp(dataset_id,'eq_xy'))
            storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/xy/xy_DataCollapse_chimq.mat';
%         end
        N_N = numel(sqrtN_vals);
        N_T = numel(T_vals);

        T_select=1:numel(T_vals);
        
        L_vals=sqrtN_vals;
        %% Bad last entries of SCF_Spin, trim that
%         for i_N = 1:N_N
%             for i_T = 1:N_T
%                 index_max=min(numel(qbin{i_N,i_T})-2, find(qbin{i_N,i_T}>.35*L_vals(i_N),1));
%                 qbin{i_N,i_T} = qbin{i_N,i_T}(1:index_max);
%                 chimq_av{i_N,i_T} = chimq_av{i_N,i_T}(1:index_max);
%             end
%         end
        
    end
    
    % Correct scaling (data not yet properly scaled from sim)
    for i_N = 1:N_N
        for i_T = 1:N_T
            q_indices=find(qbin{i_N,i_T} < pi/2);
            qbin{i_N,i_T} = qbin{i_N,i_T}(q_indices);
            chimxq_av{i_N,i_T} = chimxq_av{i_N,i_T}(q_indices) / sqrtN_vals(i_N)^2;
            chimyq_av{i_N,i_T} = chimyq_av{i_N,i_T}(q_indices) / sqrtN_vals(i_N)^2;
            chimparq_av{i_N,i_T} = chimparq_av{i_N,i_T}(q_indices) / sqrtN_vals(i_N)^2;
            chimperpq_av{i_N,i_T} = chimperpq_av{i_N,i_T}(q_indices) / sqrtN_vals(i_N)^2;
        end
    end
    eta_diff = .005;
    eta_max = .5;

%     [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);


    spline_curve_par=cell(N_N,N_T);
    spline_curve_perp=cell(N_N,N_T);
    spline_curve_xy=cell(N_N,N_T);
    spline_goodness_par=cell(N_N,N_T);
    spline_goodness_perp=cell(N_N,N_T);
    spline_goodness_xy=cell(N_N,N_T);
    spline_output_par=cell(N_N,N_T);
    spline_output_perp=cell(N_N,N_T);
    spline_output_xy=cell(N_N,N_T);
    spline_fiterr_par=zeros(N_N,N_T);
    spline_fiterr_perp=zeros(N_N,N_T);
    spline_fiterr_xy=zeros(N_N,N_T);
    
    
    for i_N = 1:length(sqrtN_vals)
        sqrtN_cur = sqrtN_vals(i_N);
        fprintf('++ Fitting. sqrtN = %d\n', sqrtN_cur);
        for i_T = 1:length(T_vals)
            T_cur = T_vals(i_T);
            q_cur = qbin{i_N,i_T};
            
            chim_cur = chimparq_av{i_N,i_T};
            [spline_curve_par{i_N,i_T}, spline_goodness_par{i_N,i_T}, spline_output_par{i_N,i_T}] = fit(q_cur(:),chim_cur(:),'smoothingspline');
            spline_fiterr_par(i_N,i_T) = sqrt(sum((chim_cur(:) - spline_curve_par{i_N,i_T}(q_cur)).^2/numel(q_cur)));
            
            chim_cur = chimperpq_av{i_N,i_T};
            [spline_curve_perp{i_N,i_T}, spline_goodness_perp{i_N,i_T}, spline_output_perp{i_N,i_T}] = fit(q_cur(:),chim_cur(:),'smoothingspline');
            spline_fiterr_perp(i_N,i_T) = sqrt(sum((chim_cur(:) - spline_curve_perp{i_N,i_T}(q_cur)).^2/numel(q_cur)));
            
            chim_cur = chimxq_av{i_N,i_T} + chimyq_av{i_N,i_T};
            [spline_curve_xy{i_N,i_T}, spline_goodness_xy{i_N,i_T}, spline_output_xy{i_N,i_T}] = fit(q_cur(:),chim_cur(:),'smoothingspline');
            spline_fiterr_xy(i_N,i_T) = sqrt(sum((chim_cur(:) - spline_curve_xy{i_N,i_T}(q_cur)).^2/numel(q_cur)));
        end
    end
    
    
    eta_testvals = 0:eta_diff:eta_max;
    eta_par_vals = 0 * T_vals;
    eta_par_err = 0 * T_vals;
    eta_perp_vals = 0 * T_vals;
    eta_perp_err = 0 * T_vals;
    eta_xy_vals = 0 * T_vals;
    eta_xy_err = 0 * T_vals;
    error_par_vals = zeros(numel(T_vals),numel(eta_testvals));
    error_perp_vals = zeros(numel(T_vals),numel(eta_testvals));
    error_xy_vals = zeros(numel(T_vals),numel(eta_testvals));
    
    % Starting index. Set to two to ignore the smallest system, which leads
    % to errors.
    i_start = 3;
    for i_T = 1:length(T_vals)
        T_cur = T_vals(i_T);
        qbins=qbin(i_start:end,i_T);
        L_cur=L_vals(i_start:end);
        
        curves=spline_curve_par(i_start:end,i_T);
        for i_eta = 1:numel(eta_testvals)
            eta=eta_testvals(i_eta);
            error_par_vals(i_T,i_eta) = CalculateCollapseError(curves,qbins,L_cur,-1,-(2-eta),0);
        end
        [errmin,err_argmin] = min(error_par_vals(i_T,:));
        eta_par_vals(i_T) = eta_testvals(err_argmin);
        if (err_argmin ~= 1 && err_argmin ~= numel(eta_testvals))
            eta_par_err(i_T) = ((error_par_vals(i_T,err_argmin+1) + error_par_vals(i_T,err_argmin-1) - 2*error_par_vals(i_T,err_argmin) )...
                / (eta_diff^2))^(-1);
        end
        
        curves=spline_curve_perp(i_start:end,i_T);
        for i_eta = 1:numel(eta_testvals)
            eta=eta_testvals(i_eta);
            error_perp_vals(i_T,i_eta) = CalculateCollapseError(curves,qbins,L_cur,-1,-(2-eta),0);
        end
        [errmin,err_argmin] = min(error_perp_vals(i_T,:));
        eta_perp_vals(i_T) = eta_testvals(err_argmin);
        if (err_argmin ~= 1 && err_argmin ~= numel(eta_testvals))
            eta_perp_err(i_T) = ((error_perp_vals(i_T,err_argmin+1) + error_perp_vals(i_T,err_argmin-1) - 2*error_perp_vals(i_T,err_argmin) )...
                / (eta_diff^2))^(-1);
        end
        
        curves=spline_curve_xy(i_start:end,i_T);
        for i_eta = 1:numel(eta_testvals)
            eta=eta_testvals(i_eta);
            error_xy_vals(i_T,i_eta) = CalculateCollapseError(curves,qbins,L_cur,-1,-(2-eta),0);
        end
        [errmin,err_argmin] = min(error_xy_vals(i_T,:));
        eta_xy_vals(i_T) = eta_testvals(err_argmin);
        if (err_argmin ~= 1 && err_argmin ~= numel(eta_testvals))
            eta_xy_err(i_T) = ((error_xy_vals(i_T,err_argmin+1) + error_xy_vals(i_T,err_argmin-1) - 2*error_xy_vals(i_T,err_argmin) )...
                / (eta_diff^2))^(-1);
        end
        
    end
   
    disp(sprintf('Saving to %s',storefilename));
    save(storefilename, 'eta_testvals',...
        'spline_curve_par','spline_fiterr_par','spline_goodness_par',...
        'eta_par_vals','eta_par_err','error_par_vals',...
        'spline_curve_perp','spline_fiterr_perp','spline_goodness_perp',...
        'eta_perp_vals','eta_perp_err','error_perp_vals',...
        'spline_curve_xy','spline_fiterr_xy','spline_goodness_xy',...
        'eta_xy_vals','eta_xy_err','error_xy_vals',...
        'sqrtN_vals','T_vals');
    
end