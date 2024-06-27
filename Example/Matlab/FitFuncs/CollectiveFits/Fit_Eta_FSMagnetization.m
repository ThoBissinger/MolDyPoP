clear all
close all
addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/FitFuncs/
addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes

%% Variable setup.




model='xy_s'; % 'xy', 'mxy' 'fmxy' 'xy_s'
rho=3.00; % 3.55, 3.00, 2.84 (2.81)
coord_num_eff = 6; % effective coordination number. Density-dependent, must be guessed.
% extractfile=sprintf('%s_reduced',integ);
if (strcmp(model,'xy'))
%     xydata=load('xy/lf0_eq.mat');
    data=load('xy/xy_dynamics_LinearTime.mat');
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
    
    
    
    runmax=250;
    storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/xy/xy_etafit_FSMag.mat';

    r_max_factor = 1; % To ignore the finite-size effect due to the periodic boundary boundary

elseif (strcmp(model,'xy_s'))
%     xydata=load('xy/lf0_eq.mat');
    data=load('xy/eq_xy_s.mat');
%     basedir="/data/scc/thobi/201207_equilibration/xy";
%     sqrtN_dirs=["anneal/sqrtN_16" "anneal/sqrtN_32", "anneal/sqrtN_64", "anneal/sqrtN_128" "scale/sqrtN_256"];
%     sampfilename="samp_eq.mat";
    basedir="/data/scc/thobi/201207_equilibration/xy_s";
    sqrtN_dirs=["anneal/sqrtN_16" "scale/sqrtN_32", "scale/sqrtN_64", "scale/sqrtN_128" "scale/sqrtN_256"];
    sampfilename="samp_eq.mat";
    
    T_dirs_cell={"T_.10" "T_.20" "T_.30" "T_.40" "T_.50" "T_.60" "T_.70" ...
        "T_.80" "T_.85" "T_.87" "T_.89" "T_.90" "T_.91" "T_.93" "T_.95" ...
        "T_.97" "T_1.00" "T_1.03" "T_1.06" "T_1.09" "T_1.10" ...
        "T_1.12" "T_1.15" "T_1.18" "T_1.20" "T_1.21" "T_1.24" "T_1.30" ......
        "T_1.40" "T_1.50" "T_1.60" "T_1.70" "T_1.80" "T_1.90" "T_2.00"};
    
    
    
    runmax=125;
    storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/xy/xy_s_etafit_FSMag.mat';

    r_max_factor = 1; % To ignore the finite-size effect due to the periodic boundary boundary
elseif (strcmp(model,'mxy'))
    data=load('mxy/rho_3.00_dynamics_LinearTime.mat');
    basedir="/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00";
    sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128" "sqrtN_256"];
%     sampfilename="samp_Dynamics_FFT.mat";
    sampfilename="samp_Dynamics.mat";
    T_dirs_cell={"T_.01", "T_.03", "T_.05", "T_.07", "T_.09", "T_.11", "T_.13", ...
        "T_.14", "T_.15", "T_.155", "T_.16", "T_.165", "T_.17", "T_.175", ...
        "T_.18", "T_.185", "T_.19", "T_.195", "T_.20", "T_.205", ...
        "T_.21", "T_.22", "T_.23", "T_.24", "T_.25", "T_.27", "T_.29", ...
        "T_.31",  "T_.33", "T_.35", "T_.37", "T_.40", "T_.43", "T_.46", "T_.49", "T_.52"};
    
    
    
    runmax=500;
    storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_etafit_FSMag.mat';
    
    r_max_factor = .5; % To ignore the finite-size effect due to the periodic boundary boundary

elseif (strcmp(model,'fmxy'))
    data=load('fmxy/fmxy_dynamics_LinearTime.mat');
    basedir="/data/scc/thobi/210715_LinearTimeSampling/fmxy";
    sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128" "sqrtN_256"];
%     sampfilename="samp_Dynamics_FFT.mat";
    sampfilename="samp_Dynamics.mat";
    T_dirs_cell={"T_.01", "T_.03", "T_.05", "T_.07", "T_.09", "T_.11", "T_.13", ...
        "T_.14", "T_.15", "T_.155", "T_.16", "T_.165", "T_.17", "T_.175", ...
        "T_.18", "T_.185", "T_.19", "T_.195", "T_.20", "T_.205", ...
        "T_.21", "T_.22", "T_.23", "T_.24", "T_.25", "T_.27", "T_.29", ...
        "T_.31",  "T_.33", "T_.35", "T_.37", "T_.40", "T_.43", "T_.46", "T_.49", "T_.52"};
    
    
    
    runmax=500;
    storefilename='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/fmxy/fmxy_etafit_FSMag.mat';
    
    r_max_factor = .5; % To ignore the finite-size effect due to the periodic boundary boundary
end
sqrtN_vals = data.('sqrtN_vals');
T_vals = data.('T_vals');
absM_av = data.('absM_av');
L_vals = data.('L_vals');

T_KT = data.('T_KT');
T_C = data.('T_C');
T_star = data.('T_star');
    


N_N = length(sqrtN_vals);
N_T = length(T_vals);



%% 1 critical exponent eta and xi from SCF_r
section_ind = 1;
T_size=size(T_vals);
fitobs=cell(T_size);

eta_vals=zeros(T_size);
c_vals=zeros(T_size);

eta_ci=zeros(N_T,2);
eta_sd=zeros(T_size);

fit_confidentiality = .95;
threshold = 0; % Could be made N-dependent...
T_dirs = T_dirs_cell;
for i_T = 1:numel(T_dirs)
    T_cur = T_vals(i_T);
    disp("  T = " + T_cur);
    T_string=split(T_dirs{i_T},"_");

    M_cur=cell2mat(absM_av(:,i_T));
    fitob = fit_eta_Magnetization_FS(M_cur,L_vals);
    
    fitobs{i_T} = fitob;
    eta_vals(i_T) = fitob.eta;
    c_vals = fitob.c;
    ci=confint(fitob,fit_confidentiality);
    eta_ci(i_T,:)=ci(:,end);
    eta_sd(i_T) = sqrt(runmax * (eta_ci(i_T,2) - eta_ci(i_T,1)).^2 / (2 * 1.96));

end

disp(storefilename);
save(storefilename, 'fitobs', 'fit_confidentiality', ...
    'eta_vals', 'eta_ci','eta_sd', 'c_vals', ...
    'L_vals','T_vals','sqrtN_vals');

