clear all
close all
addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/FitFuncs/
addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes

%% Variable setup.
mxydata=load('mxy/rho_3.00_eq.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat');
xydata=load('xy/lf0_eq.mat'); xyfit=load('xy/lf0_CritExpFit.mat');


model='mxy'; % 'xy', 'mxy'
rho=3.00; % 3.55, 3.00, 2.84 (2.81)
coord_num_eff = 6; % effective coordination number. Density-dependent, must be guessed.
% extractfile=sprintf('%s_reduced',integ);
if (strcmp(model,'xy'))
    storagebase=sprintf('%s/%s',model,integ);
    suffix=sprintf('%s',integ);
    load(sprintf('%s_reduced.mat',storagebase))
    load(sprintf('%s_noqfull.mat',storagebase))
    load(sprintf('%s_FSFitting.mat',storagebase))
    load(sprintf('%s_qreduced.mat',storagebase))
elseif (strcmp(model,'mxy'))
    storagebase=sprintf('%s/rho_%1.2f',model,rho);
    suffix=sprintf('rho_%.2f',rho);

    sqrtN_vals = mxydata.('sqrtN_vals');
    T_vals = mxydata.('T_vals');
    SCF_Spin_av = mxydata.('SCF_Spin_av');
    rbin = mxydata.('rbin');
    T_KT = mxydata.('T_KT');
    T_C = mxydata.('T_C');
    T_star = mxydata.('T_star');
    
end
storefilename=sprintf('%s_CritExpFit.mat',storagebase);
plotfoldname=sprintf('%s/plots/CritExpFit/%s',model,suffix);


offset_switch=false;
% load(sprintf('%s_full',integ))


N_N = length(sqrtN_vals);
N_T = length(T_vals);

% load(storefilename)
if (strcmp(model,'xy'))
    sqrtN_total = 256;
    lattice_spacing = 2;
elseif (strcmp(model,'mxy'))
    sqrtN_total = 256;
    lattice_spacing = 1/sqrt(rho);
end
i_N_total = find(sqrtN_vals == sqrtN_total);
r_vals_total=cell2mat(rbin(i_N_total,1));


T_KT_select = find(T_vals > T_KT(i_N_total),1)-2 : find(T_vals > T_C(i_N_total),1)+1;
T_KT_select(end+1) = T_KT_select(end)+5;
% T_KT_select(end+1) = min(T_KT_select(end)+5,length(T_vals));

i_Tmin=find(T_vals<T_KT(i_N_total),1,'last');
i_Tmax=find(T_vals<T_C(i_N_total),1,'last');
% T_LessKT_select=[i_Tmin, round(.5*(i_Tmin + i_Tmax)), i_Tmax];
T_LessKT_select=1 : i_Tmax + 4;



%% 1 critical exponent eta and xi from SCF_r
section_ind = 1;
T_select=T_KT_select;
eta_SCF_initguess=zeros(length(sqrtN_vals),length(T_vals));
param_SCFSpin_Exp=cell(N_N,N_T);
param_SCFSpin_Pow=cell(N_N,N_T);

eta_vals=zeros(N_N,N_T);
xi_vals=zeros(N_N,N_T);

eta_err=zeros(N_N,N_T,2);
xi_err=zeros(N_N,N_T,2);

fit_confidentiality = .95;
threshold = 0; % Could be made N-dependent...
for i_N = 1:length(sqrtN_vals)
    sqrtN_cur = sqrtN_vals(i_N);
    for i_T = 1:length(T_vals)
        T_cur = T_vals(i_T);

        eta_SCF_initguess(i_T,i_N) = T_cur*coord_num_eff/(8*pi); % Spin wave value.

        r_vals = cell2mat(rbin(i_N,i_T));
        spin_cur=cell2mat(SCF_Spin_av(i_N, i_T));
% 
%         r_min_ind=find(r_vals > 2*lattice_spacing,1);
%         r_max_ind=find(r_vals > .8 * r_vals(end),1);
%         r_select=r_min_ind:r_max_ind;
%         r_cur = r_vals(r_select);
        r_cur = r_vals;
        r_select = 1:numel(r_cur);

        spin_cur=SCF_Spin_av{i_N, i_T}(r_select);

        if (length(r_cur) > 3)
            fitob=fit_PowSCF(r_cur,spin_cur,r_cur(1),r_cur(end),offset_switch);
            param_SCFSpin_Pow{i_N, i_T} = fitob;
            eta_vals(i_N,i_T) = fitob.eta;
            ci=confint(fitob,fit_confidentiality);
            eta_err(i_N,i_T,:)=ci(:,end);

            fitob=fit_ExpSCF(r_cur,spin_cur,r_cur(1),r_cur(end));
            param_SCFSpin_Exp{i_N, i_T} = fitob;
            xi_vals(i_N,i_T) = fitob.xi;
            ci=confint(fitob,fit_confidentiality);
            xi_err(i_N,i_T,:)=ci(:,end);
        else
            eta_vals(i_N,i_T) = 0;
            eta_err(i_N,i_T,:)= 0;
            xi_vals(i_N,i_T) = 0;
            xi_err(i_N,i_T,:)= 0;
        end


    end
end

    save(storefilename, 'param_SCFSpin_Exp', 'param_SCFSpin_Pow', ...
        'eta_vals', 'xi_vals','eta_err', 'xi_err', 'fit_confidentiality',...
        'eta_SCF_initguess','T_vals','sqrtN_vals');

