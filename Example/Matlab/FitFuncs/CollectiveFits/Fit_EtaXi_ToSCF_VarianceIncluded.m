clear all
close all
addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/FitFuncs/
addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes

%% Variable setup.




model='fmxy'; % 'xy', 'mxy' 'fmxy'
rho=3.00; % 3.55, 3.00, 2.84 (2.81)
coord_num_eff = 6; % effective coordination number. Density-dependent, must be guessed.
% extractfile=sprintf('%s_reduced',integ);
if (strcmp(model,'xy'))
%     xydata=load('xy/lf0_eq.mat');
    xydata=load('xy/xy_dynamics_LinearTime.mat');
%     basedir="/data/scc/thobi/201207_equilibration/xy";
%     sqrtN_dirs=["anneal/sqrtN_16" "anneal/sqrtN_32", "anneal/sqrtN_64", "anneal/sqrtN_128" "scale/sqrtN_256"];
%     sampfilename="samp_eq.mat";
    basedir="/data/scc/thobi/210715_LinearTimeSampling/xy";
    sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128" "sqrtN_256"];
    sampfilename="samp_Dynamics.mat";
    
    T_dirs_cell={"T_.10" "T_.20" "T_.30" "T_.40" "T_.50" "T_.60" "T_.70" "T_.80" "T_.90" "T_1.00" "T_1.10" "T_1.20" "T_1.30"...
        "T_1.40" "T_1.42" "T_1.44" "T_1.46" "T_1.48" "T_1.50" "T_1.52" "T_1.54" "T_1.56" "T_1.58" "T_1.60" "T_1.62" "T_1.64" "T_1.66"...
        "T_1.70" "T_1.75" "T_1.80" "T_1.85" "T_1.90" "T_1.95" "T_2.00" "T_2.10" "T_2.20" "T_2.30" "T_2.40" "T_2.50" "T_2.60" "T_2.70"...
        "T_2.80" "T_2.90" "T_3.00"};
    
    
    sqrtN_vals = xydata.('sqrtN_vals');
    T_vals = xydata.('T_vals');
    
    T_KT = xydata.('T_KT');
    T_C = xydata.('T_C');
    T_star = xydata.('T_star');
    
    runmax=250;
    storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/xy/xy_CritExpFit.mat';

    r_max_factor = 1; % To ignore the finite-size effect due to the periodic boundary boundary
elseif (strcmp(model,'mxy'))
    mxydata=load('mxy/rho_3.00_dynamics_LinearTime.mat');
    basedir="/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00";
    sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128" "sqrtN_256"];
%     sampfilename="samp_Dynamics_FFT.mat";
    sampfilename="samp_Dynamics.mat";
    T_dirs_cell={"T_.01", "T_.03", "T_.05", "T_.07", "T_.09", "T_.11", "T_.13", ...
        "T_.14", "T_.15", "T_.155", "T_.16", "T_.165", "T_.17", "T_.175", ...
        "T_.18", "T_.185", "T_.19", "T_.195", "T_.20", "T_.205", ...
        "T_.21", "T_.22", "T_.23", "T_.24", "T_.25", "T_.27", "T_.29", ...
        "T_.31",  "T_.33", "T_.35", "T_.37", "T_.40", "T_.43", "T_.46", "T_.49", "T_.52"};
    
    
    sqrtN_vals = mxydata.('sqrtN_vals');
    T_vals = mxydata.('T_vals');
    
    T_KT = mxydata.('T_KT');
    T_C = mxydata.('T_C');
    T_star = mxydata.('T_star');
    
    runmax=500;
    storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_CritExpFit.mat';
    
    r_max_factor = .5; % To ignore the finite-size effect due to the periodic boundary boundary

elseif (strcmp(model,'fmxy'))
    fmxydata=load('fmxy/fmxy_dynamics_LinearTime.mat');
    basedir="/data/scc/thobi/210715_LinearTimeSampling/fmxy";
    sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128" "sqrtN_256"];
%     sampfilename="samp_Dynamics_FFT.mat";
    sampfilename="samp_Dynamics.mat";
    T_dirs_cell={"T_.01", "T_.03", "T_.05", "T_.07", "T_.09", "T_.11", "T_.13", ...
        "T_.14", "T_.15", "T_.155", "T_.16", "T_.165", "T_.17", "T_.175", ...
        "T_.18", "T_.185", "T_.19", "T_.195", "T_.20", "T_.205", ...
        "T_.21", "T_.22", "T_.23", "T_.24", "T_.25", "T_.27", "T_.29", ...
        "T_.31",  "T_.33", "T_.35", "T_.37", "T_.40", "T_.43", "T_.46", "T_.49", "T_.52"};
    
    
    sqrtN_vals = fmxydata.('sqrtN_vals');
    T_vals = fmxydata.('T_vals');
    
    T_KT = fmxydata.('T_KT');
    T_C = fmxydata.('T_C');
    T_star = fmxydata.('T_star');
    
    runmax=500;
    storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/fmxy/fmxy_CritExpFit.mat';
    
    r_max_factor = .5; % To ignore the finite-size effect due to the periodic boundary boundary
end

offset_switch=false;

N_N = length(sqrtN_vals);
N_T = length(T_vals);



%% 1 critical exponent eta and xi from SCF_r
section_ind = 1;
% T_select=T_KT_select;
eta_SCF_initguess=zeros(length(sqrtN_vals),length(T_vals));
param_SCFSpin_Exp=cell(N_N,N_T);
param_SCFSpin_Pow=cell(N_N,N_T);

eta_vals=zeros(N_N,N_T);
xi_vals=zeros(N_N,N_T);

eta_ci=zeros(N_N,N_T,2);
eta_sd=zeros(N_N,N_T);
xi_ci=zeros(N_N,N_T,2);
xi_sd=zeros(N_N,N_T);

fit_confidentiality = .95;
threshold = 0; % Could be made N-dependent...
for i_N = 1:numel(sqrtN_dirs)
    sqrtN_cur = sqrtN_vals(i_N);
    sqrtN_string=split(sqrtN_dirs(i_N),"_");
    T_dirs = T_dirs_cell;
    for i_T = 1:numel(T_dirs)
        T_cur = T_vals(i_T);
        disp("sqrtN = " + sqrtN_cur + "  T = " + T_cur);
        T_string=split(T_dirs{i_T},"_");
        cur_collectfile=sprintf('%s/%s/%s/%s',basedir,sqrtN_dirs(i_N),T_dirs{i_T},sampfilename);
        load(cur_collectfile,'rbin_collect','SCF_Spin_av_collect');
        if ( numel(SCF_Spin_av_collect(:,1)) ~= runmax ) % Warns in length mismatch
            SCF_Spin_av_collect = SCF_Spin_av_collect(1:runmax,:);
            rbin_collect = rbin_collect(1:runmax,:);
            fprintf('Length problem in %s\n',cur_collectfile)
        end
%         curpath=sprintf('%s/%s/%s',basedir,sqrtN_dirs(i_N),T_dirs{i_T});
%                 pathname = sprintf('%s/sqrtN_%d/T_%s', pathbase, sqrtN,T_strings{iT});

%         fullfilename = sprintf('%s/%s',curpath,collectfilename);
%         disp(fullfilename)

        eta_SCF_initguess(i_T,i_N) = T_cur*coord_num_eff/(8*pi); % Spin wave value.

%         r_vals = rbin_collect;
%         spin_cur=SCF_Spin_av_collect;
        nonzero_ind=find(SCF_Spin_av_collect);
        
        r_vals = rbin_collect(nonzero_ind);
        spin_cur=SCF_Spin_av_collect(nonzero_ind);
        
%         if (length(r_cur) > 3)
            fitob=fit_PowSCF(r_vals,spin_cur,min(r_vals),r_max_factor * max(r_vals),offset_switch);
            param_SCFSpin_Pow{i_N, i_T} = fitob;
            eta_vals(i_N,i_T) = fitob.eta;
            ci=confint(fitob,fit_confidentiality);
            eta_ci(i_N,i_T,:)=ci(:,end);
            eta_sd(i_N,i_T) = sqrt(runmax * (eta_ci(i_N,i_T,2) - eta_ci(i_N,i_T,1)).^2 / (2 * 1.96));

            fitob=fit_ExpSCF(r_vals,spin_cur,min(r_vals),r_max_factor * max(r_vals));
            param_SCFSpin_Exp{i_N, i_T} = fitob;
            xi_vals(i_N,i_T) = fitob.xi;
            ci=confint(fitob,fit_confidentiality);
            xi_ci(i_N,i_T,:)=ci(:,end);
            xi_sd(i_N,i_T) = sqrt(runmax * (xi_ci(i_N,i_T,2) - xi_ci(i_N,i_T,1)).^2 / (2 * 1.96));


    end
end

    save(storefilename, 'param_SCFSpin_Exp', 'param_SCFSpin_Pow', ...
        'eta_vals', 'xi_vals','eta_ci','eta_sd', 'xi_ci', 'xi_sd', 'fit_confidentiality',...
        'eta_SCF_initguess','T_vals','sqrtN_vals');

