init%% Initializing

clear all
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;

% plot_type=["w", "mperp", "te","mpar"] %% , "m"];
corr_fct=["mperp"];
% set_identifier = 'LinearTime';
set_identifier = 'AdjustedTime';

for i_model = [4]

    % figure(i_model)
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        
        if (strcmp(set_identifier,'LinearTime'))
            data=matfile('mxy/rho_3.00_dynamics_LinearTime.mat');
        else
            data=matfile('mxy/mxy_dynamics_AdjustedTime_SmallN.mat');
        end
        
        storefilebase='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_TCF_q_fit';
        storefilesuffix=set_identifier;
        

        L_vals=[9.25,18.5,37,74,148];
    elseif (i_model == 2)
        curmodel="xy";
        curtitle="SXY model";
        
        data=matfile('xy/xy_dynamics_LinearTime.mat');
        
%         storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/xy/xy_TCF_q_fit_LinearTime.mat';
        storefilebase='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/xy/xy_TCF_q_fit';
        storefilesuffix='LinearTime';
%         storefilename=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/xy/xy_TCF_q_fit_%s_LinearTime.mat',corr_fct);
        
        L_vals=[16,32,64,128,256];
    elseif (i_model == 3)
        curmodel="fmxy";
        curtitle="FMXY model";
        
        
        if (strcmp(set_identifier,'LinearTime'))
            data=matfile('fmxy/fmxy_dynamics_LinearTime.mat');
        else
            data=matfile('fmxy/fmxy_dynamics_AdjustedTime.mat');
        end
        
%         storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/fmxy/fmxy_TCF_q_fit_LinearTime.mat';
%         storefilename=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/fmxy/fmxy_TCF_q_fit_%s_LinearTime.mat',corr_fct);
        storefilebase='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/fmxy/fmxy_TCF_q_fit';
        storefilesuffix=set_identifier;
        
        
        L_vals=[9.25,18.5,37,74,148];
     elseif (i_model == 4)
        curmodel="xy_s";
        data=matfile('xy/xy_s_dynamics_AdjustedTime.mat');
        
%         storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/fmxy/fmxy_TCF_q_fit_LinearTime.mat';
%         storefilename=sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/fmxy/fmxy_TCF_q_fit_%s_LinearTime.mat',corr_fct);
        storefilebase='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/xy/xy_s_TCF_q_fit';
        storefilesuffix=set_identifier;
        
        L_vals=[16, 32, 64, 128, 256];
    end
    sqrtN_vals=data.('sqrtN_vals');
    T_vals=data.('T_vals');
    N_N = numel(sqrtN_vals);
%     gmperpmperp=data.('gmperpmperp');
    ACF_q0_M=data.('ACF_q0_M');
    absM_av=data.('absM_av');
    TCF_times=data.('TCF_times');
    
    N_T = numel(T_vals);
    for i_cfct = 1:numel(corr_fct)
        corr_fct_cur = corr_fct(i_cfct);
        qbin = data.('qbin');
        storefilename=sprintf('%s_%s_%s.mat',storefilebase,corr_fct_cur,storefilesuffix);
%     if (plot_type_cur == "m")
%         cf = gxx{i_N,i_T}(q_indices) + gyy{i_N,i_T}(q_indices);
%     else
        if (corr_fct_cur == "mpar")
            cf = data.('gmparmpar');
        elseif (corr_fct_cur == "mperp")
            cf = data.('gmperpmperp');
        elseif (corr_fct_cur == "te")
            cf = data.('gtt');
        elseif (corr_fct_cur == "w")
            cf = data.('gww');
        end
    
        n_period = 6;
        weightexp = .25;
    

    %% Fitting

    %     [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
        qbin_with_zero_clipped=cell(N_N,N_T);
        
        param_TCFSpin_omega_1_q_DO=cell(N_N,N_T);
        fitobs_TCFSpin_omega_1_q_DO=cell(N_N,N_T);
        fitobs_gamma_with_offset=cell(N_N,N_T);
        fitobs_gamma=cell(N_N,N_T);
        fitobs_omega_1=cell(N_N,N_T);
        
        param_TCFSpin_omega_0_q_DO=cell(N_N,N_T);
    
        gamma_a=zeros(N_N,N_T);
        gamma_sigma=zeros(N_N,N_T);
        gamma_a_with_offset=zeros(N_N,N_T);
        gamma_sigma_with_offset=zeros(N_N,N_T);
        gamma_c_with_offset=zeros(N_N,N_T);
        omega_1_a=zeros(N_N,N_T);
        omega_1_sigma=zeros(N_N,N_T);
        
    %     eta_vals=zeros(N_N,N_T);
    %     xi_vals=zeros(N_N,N_T);
    % 
    %     eta_err=zeros(N_N,N_T,2);
    %     xi_err=zeros(N_N,N_T,2);
    
        fit_confidentiality = .95;
        threshold = 0; % Could be made N-dependent...
        for i_N = 1:length(sqrtN_vals)
            sqrtN_cur = sqrtN_vals(i_N);
            
            fprintf('++ Fitting. sqrtN = %d\n', sqrtN_cur);
            for i_T = 1:length(T_vals)
                T_cur = T_vals(i_T);
                fprintf('+++ sqrtN = %d, T = %.3f,', sqrtN_cur,T_cur);
                n_q_old = numel(qbin{i_N,i_T}); % is the same as n_q currently, but just to be sure
                q_vals = [0,qbin{i_N,i_T}(1:end-1)]; % Removes the last entry, it is always an empty bin (should be changed in code)
                n_q = numel(q_vals);
                qbin_with_zero_clipped{i_N,i_T} = q_vals;
                
                param_TCFSpin_omega_1_q_DO{i_N,i_T}=zeros(1,2*length(q_vals));
                coeffs_omega_1_cur=zeros(1,2*length(q_vals));
                coeffs_omega_0_cur=zeros(1,2*length(q_vals));
                cf_cur = cf{i_N,i_T};
                t_vals=TCF_times{i_N,i_T};
                fitobs_cur = cell(length(q_vals));
                for i_q = 1:length(q_vals) 
                    q_cur = q_vals(i_q);
                    fprintf(' q = %.3f, ', q_cur);
                    if q_cur == 0
                        TCF_Spin_cur = ACF_q0_M{i_N,i_T};
                    else
                        TCF_Spin_cur = cf_cur(i_q-1:n_q_old:end);
                        % Note that this ignores the last entry
                    end
                    
                    % For better agreement (mostly towards the rear) maybe smoothen
    %                 TCF_Spin_cur=smoothen(TCF_Spin_cur,3,100);
                    coeffs_omega_1_cur(2*(i_q-1)+(1:2))=fit_DampedOscillator_RealSpace(t_vals,real(TCF_Spin_cur),n_period,weightexp,'omega_1');
                    coeffs_omega_0_cur(2*(i_q-1)+(1:2))=fit_DampedOscillator_RealSpace(t_vals,real(TCF_Spin_cur),n_period,weightexp,'omega_0');
    %!     ALTERNATIVE FIT. SO FAR BAD RESULTS, BUT MAYBE STILL IMPROVABLE.
    %                 fitobs_cur{i_q}=fit_DampedOscillator_alternative(t_vals,real(TCF_Spin_cur));
                end
                fprintf('\n');
                
                param_TCFSpin_omega_1_q_DO{i_N,i_T}=coeffs_omega_1_cur;
                param_TCFSpin_omega_0_q_DO{i_N,i_T}=coeffs_omega_0_cur;
    %             fitobs_TCFSpin_q_DO{i_N,i_T}=fitobs_cur;
    
                % !! Fitting of coefficients to power laws.
                % ind_end: We must ignore the last value, it is corrupt 
                % (has to do with the binning). This is ok by the 
                % redefiniton of q_vals above. Later values may also 
                % contain binned quantities, leading to worse agreement.
                ind_end = min(max(ceil(length(q_vals)/2),8),length(q_vals));
                % ind_start: first index is the q=0 value. This will not be
                % included in the fit, it is of a different nature...
                ind_start = 2; 
                
                f_FIT = fittype('a*(x)^(sigma)+c');
                curcoeff=coeffs_omega_1_cur(1:2:end);
                curfitob = fit(q_vals(ind_start:ind_end)', curcoeff(ind_start:ind_end)',f_FIT,...
                    'StartPoint',[curcoeff(ind_start)/q_vals(ind_start)^2,2,0],'lower',[0,0,0]); % gamma fit with offset
                gamma_a_with_offset(i_N,i_T) = curfitob.a;
                gamma_sigma_with_offset(i_N,i_T) = curfitob.sigma;
                gamma_c_with_offset(i_N,i_T) = curfitob.c;
                fitobs_gamma_with_offset{i_N,i_T} = curfitob;
    
                f_FIT = fittype('a*(x)^(sigma)');
                curcoeff=coeffs_omega_1_cur(1:2:end);
                curfitob = fit(q_vals(ind_start:ind_end)', curcoeff(ind_start:ind_end)',f_FIT,...
                    'StartPoint',[curcoeff(ind_start)/q_vals(ind_start)^2,2],'lower',[0,0]); % gamma fit without offset
                gamma_a(i_N,i_T) = curfitob.a;
                gamma_sigma(i_N,i_T) = curfitob.sigma;
                fitobs_gamma{i_N,i_T} = curfitob;
                
                f_FIT = fittype('a*(x)^(sigma)');
                curcoeff=coeffs_omega_1_cur(2:2:end);
                curfitob = fit(q_vals(ind_start:ind_end)', curcoeff(ind_start:ind_end)',f_FIT,...
                    'StartPoint',[curcoeff(ind_start)/q_vals(ind_start),1],'lower',[0,0]); % omega fit
                omega_1_a(i_N,i_T) = curfitob.a;
                omega_1_sigma(i_N,i_T) = curfitob.sigma;
                fitobs_omega_1{i_N,i_T} = curfitob;
                
            end
        end

        %% Printing, Cleaning, Saving
        qbin = qbin_with_zero_clipped;
        clear r_vals i_N_total sqrtN_total;

        fprintf('\n\n');
        fprintf('gamma\n N/T');
        fprintf('    %.3f', T_vals); 
        fprintf('\n');
        for i_N = 1:numel(sqrtN_vals)
            fprintf('%4d',sqrtN_vals(i_N));
            for i_T = 1:numel(T_vals)
                fprintf(' %.6f', param_TCFSpin_omega_1_q_DO{i_N,i_T}(3)); 
            end
            fprintf('\n');
        end
        fprintf('\n');

        fprintf('\n\n');
        fprintf('omega_1\n N/T');
        fprintf('    %.3f', T_vals); 
        fprintf('\n');
        for i_N = 1:numel(sqrtN_vals)
            fprintf('%4d',sqrtN_vals(i_N));
            for i_T = 1:numel(T_vals)
                fprintf(' %.6f', param_TCFSpin_omega_1_q_DO{i_N,i_T}(4)); 
            end
            fprintf('\n');
        end
        fprintf('\n');

        fprintf('\n\n');
        fprintf('1/gamma\n N/T');
        fprintf('     %.3f', T_vals); 
        fprintf('\n');
        for i_N = 1:numel(sqrtN_vals)
            fprintf('%4d',sqrtN_vals(i_N));
            for i_T = 1:numel(T_vals)
                fprintf(' %9.2f', 1/param_TCFSpin_omega_1_q_DO{i_N,i_T}(3)); 
            end
            fprintf('\n');
        end
        fprintf('\n');

        fprintf('\n\n');
        fprintf('1/omega_1\n N/T');
        fprintf('    %.3f', T_vals); 
        fprintf('\n');
        for i_N = 1:numel(sqrtN_vals)
            fprintf('%4d',sqrtN_vals(i_N));
            for i_T = 1:numel(T_vals)
                fprintf(' %8.2f', 1/param_TCFSpin_omega_1_q_DO{i_N,i_T}(4)); 
            end
            fprintf('\n');
        end
        fprintf('\n');

        disp(storefilename);
        save(storefilename, 'T_vals', 'sqrtN_vals', 'qbin', ...
            'param_TCFSpin_omega_1_q_DO','fitobs_TCFSpin_omega_1_q_DO',...
            'fitobs_gamma','gamma_a','gamma_sigma',...
            'param_TCFSpin_omega_0_q_DO',...
            'fitobs_gamma_with_offset','gamma_a_with_offset','gamma_sigma_with_offset','gamma_c_with_offset',...
            'fitobs_omega_1','omega_1_a','omega_1_sigma');
    end
end