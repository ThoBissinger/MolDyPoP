clear all;
close all;

addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;
fig_select=[1:3];
fig_base = "/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/FFT";
sim_data_dir="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis";
models = ["mxy", "xy", "fmxy"];
% models = ["mxy", "xy"];
% curmodel = "mxy";
for i_model = 1
    curmodel = models(i_model);
    % mxydata=load('mxy/rho_3.00_integ.mat');
    % mxydata=load('mxy/rho_3.00_dynamics_better_q.mat'); mxyfit=load('mxy/rho_3.00_CritExpFit.mat');
    if (curmodel == "mxy")
%         data=load('mxy/rho_3.00_dynamics_LinearTime.mat');
%         fit_FS=load('mxy/rho_3.00_DataCollapse_SCF.mat');
        collectfilepath=sprintf('%s/mxy/rho_3.00_FFT.mat',sim_data_dir);
        FFT_Collectfile='samp_Dynamics_FFT';
        data_basedir='/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00';
        T_dirs_cell={"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
        L_vals = [9.25,18.5,37,74,148];
        runmax=500;
    elseif (curmodel == "xy")
%         data=load('xy/xy_dynamics_LinearTime.mat');
%         fit_FS=load('xy/xy_DataCollapse_SCF.mat');
        collectfilepath=sprintf('%s/xy/xy_FFT.mat',sim_data_dir);
        FFT_Collectfile='samp_Dynamics_FFT';
        data_basedir='/data/scc/thobi/210715_LinearTimeSampling/xy';
        T_dirs_cell={"T_.10", "T_.20", "T_.30", "T_.40", "T_.50", "T_.60", "T_.70", "T_.80", "T_.90", "T_1.00"  "T_1.10" "T_1.20" "T_1.30" "T_1.40" "T_1.42" "T_1.44" "T_1.46" "T_1.48" "T_1.50" "T_1.52" "T_1.54" "T_1.56" "T_1.58" "T_1.60" "T_1.62" "T_1.64" "T_1.66" "T_1.70" "T_1.75" "T_1.80" "T_1.85" "T_1.90" "T_1.95" "T_2.00" "T_2.10" "T_2.20" "T_2.30" "T_2.40" "T_2.50" "T_2.60" "T_2.70" "T_2.80" "T_2.90" "T_3.00"};
        L_vals=[16,32,64,128,256];
        runmax=250;
    elseif (curmodel == "fmxy")
%         data=load('fmxy/fmxy_dynamics_LinearTime.mat');
%         fit_FS=load('fmxy/fmxy_DataCollapse_SCF.mat');
        collectfilepath=sprintf('%s/fmxy/fmxy_FFT.mat',sim_data_dir);
        FFT_Collectfile='samp_Dynamics_FFT';
        data_basedir='/data/scc/thobi/210715_LinearTimeSampling/fmxy';
        T_dirs_cell={"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
        L_vals = [9.25,18.5,37,74,148];
        runmax=500;
    end
    sqrtN_vals=[16,32,64,128,256];
    N_N=numel(sqrtN_vals);
    N_T=numel(T_dirs_cell);


    Sxx_noiselevel = cell(N_N,N_T);
    Sxx_h_0 = cell(N_N,N_T);
    Sxx_omega_0 = cell(N_N,N_T);
    Sxx_gamma = cell(N_N,N_T);
    Sxx_fitob = cell(N_N,N_T);
    Syy_noiselevel = cell(N_N,N_T);
    Syy_h_0 = cell(N_N,N_T);
    Syy_omega_0 = cell(N_N,N_T);
    Syy_gamma = cell(N_N,N_T);
    Syy_fitob = cell(N_N,N_T);
    Smperpmperp_noiselevel = cell(N_N,N_T);
    Smperpmperp_h_0 = cell(N_N,N_T);
    Smperpmperp_omega_0 = cell(N_N,N_T);
    Smperpmperp_gamma = cell(N_N,N_T);
    Smperpmperp_fitob = cell(N_N,N_T);
    Smparmpar_noiselevel = cell(N_N,N_T);
    Smparmpar_h_0 = cell(N_N,N_T);
    Smparmpar_omega_0 = cell(N_N,N_T);
    Smparmpar_gamma = cell(N_N,N_T);
    Smparmpar_fitob = cell(N_N,N_T);
    Stt_noiselevel = cell(N_N,N_T);
    Stt_h_0 = cell(N_N,N_T);
    Stt_omega_0 = cell(N_N,N_T);
    Stt_gamma = cell(N_N,N_T);
    Stt_fitob = cell(N_N,N_T);
    Sww_noiselevel = cell(N_N,N_T);
    Sww_h_0 = cell(N_N,N_T);
    Sww_omega_0 = cell(N_N,N_T);
    Sww_gamma = cell(N_N,N_T);
    Sww_fitob = cell(N_N,N_T);

%     Sxx=cell(N_N,N_T);
%     Syy=cell(N_N,N_T);
%     Smparmpar=cell(N_N,N_T);
%     Smperpmperp=cell(N_N,N_T);
%     Stt=cell(N_N,N_T);
%     Sww=cell(N_N,N_T);
% 
%     Sxx_noiselevel = cell(N_N,N_T);
%     Sxx_h_0 = cell(N_N,N_T);
%     Sxx_omega_0 = cell(N_N,N_T);
%     Sxx_gamma = cell(N_N,N_T);
%     Sxx_fitob = cell(N_N,N_T);
% 
%     Syy_noiselevel = cell(N_N,N_T);
%     Syy_h_0 = cell(N_N,N_T);
%     Syy_omega_0 = cell(N_N,N_T);
%     Syy_gamma = cell(N_N,N_T);
%     Syy_fitob = cell(N_N,N_T);
% 
%     Smperpmperp_noiselevel = cell(N_N,N_T);
%     Smperpmperp_h_0 = cell(N_N,N_T);
%     Smperpmperp_omega_0 = cell(N_N,N_T);
%     Smperpmperp_gamma = cell(N_N,N_T);
%     Smperpmperp_fitob = cell(N_N,N_T);
% 
%     Smparmpar_noiselevel = cell(N_N,N_T);
%     Smparmpar_h_0 = cell(N_N,N_T);
%     Smparmpar_omega_0 = cell(N_N,N_T);
%     Smparmpar_gamma = cell(N_N,N_T);
%     Smparmpar_fitob = cell(N_N,N_T);
% 
%     Stt_noiselevel = cell(N_N,N_T);
%     Stt_h_0 = cell(N_N,N_T);
%     Stt_omega_0 = cell(N_N,N_T);
%     Stt_gamma = cell(N_N,N_T);
%     Stt_fitob = cell(N_N,N_T);
% 
%     Sww_noiselevel = cell(N_N,N_T);
%     Sww_h_0 = cell(N_N,N_T);
%     Sww_omega_0 = cell(N_N,N_T);
%     Sww_gamma = cell(N_N,N_T);
%     Sww_fitob = cell(N_N,N_T);
% 
%     qbin=cell(N_N,N_T);
%     om_vals=cell(N_N,N_T);
%     averaging_times=cell(N_N,N_T);
    
    
    for i_N = 1:N_N
        sqrtN=sqrtN_vals(i_N);
        for i_T = 1:N_T
            fprintf('%d, %d\n',i_N,i_T);
            cur_collectfile=sprintf('%s/sqrtN_%d/%s/%s',data_basedir,sqrtN,...
                T_dirs_cell{i_T},FFT_Collectfile);
            curmat=matfile(cur_collectfile);
            
            Sxx_noiselevel{i_N,i_T} = curmat.Sxx_noiselevel;
            Sxx_h_0{i_N,i_T} = curmat.Sxx_h_0;
            Sxx_omega_0{i_N,i_T} = curmat.Sxx_omega_0;
            Sxx_gamma{i_N,i_T} = curmat.Sxx_gamma;
%             Sxx_fitob{i_N,i_T} = curmat.Sxx_fitob;
            Syy_noiselevel{i_N,i_T} = curmat.Syy_noiselevel;
            Syy_h_0{i_N,i_T} = curmat.Syy_h_0;
            Syy_omega_0{i_N,i_T} = curmat.Syy_omega_0;
            Syy_gamma{i_N,i_T} = curmat.Syy_gamma;
%             Syy_fitob{i_N,i_T} = curmat.Syy_fitob;
            Smperpmperp_noiselevel{i_N,i_T} = curmat.Smperpmperp_noiselevel;
            Smperpmperp_h_0{i_N,i_T} = curmat.Smperpmperp_h_0;
            Smperpmperp_omega_0{i_N,i_T} = curmat.Smperpmperp_omega_0;
            Smperpmperp_gamma{i_N,i_T} = curmat.Smperpmperp_gamma;
%             Smperpmperp_fitob{i_N,i_T} = curmat.Smperpmperp_fitob;
            Smparmpar_noiselevel{i_N,i_T} = curmat.Smparmpar_noiselevel;
            Smparmpar_h_0{i_N,i_T} = curmat.Smparmpar_h_0;
            Smparmpar_omega_0{i_N,i_T} = curmat.Smparmpar_omega_0;
            Smparmpar_gamma{i_N,i_T} = curmat.Smparmpar_gamma;
%             Smparmpar_fitob{i_N,i_T} = curmat.Smparmpar_fitob;
            Stt_noiselevel{i_N,i_T} = curmat.Stt_noiselevel;
            Stt_h_0{i_N,i_T} = curmat.Stt_h_0;
            Stt_omega_0{i_N,i_T} = curmat.Stt_omega_0;
            Stt_gamma{i_N,i_T} = curmat.Stt_gamma;
%             Stt_fitob{i_N,i_T} = curmat.Stt_fitob;
            Sww_noiselevel{i_N,i_T} = curmat.Sww_noiselevel;
            Sww_h_0{i_N,i_T} = curmat.Sww_h_0;
            Sww_omega_0{i_N,i_T} = curmat.Sww_omega_0;
            Sww_gamma{i_N,i_T} = curmat.Sww_gamma;
%             Sww_fitob{i_N,i_T} = curmat.Sww_fitob;

            

%             Sxx{i_N,i_T} = curmat.Sxx;
%             Syy{i_N,i_T} = curmat.Syy;
%             Smperpmperp{i_N,i_T} = curmat.Smperpmperp;
%             Smparmpar{i_N,i_T} = curmat.Smparmpar;
%             Stt{i_N,i_T} = curmat.Stt;
%             Sww{i_N,i_T} = curmat.Sww;
% 
%             Sxx_noiselevel{i_N,i_T} = curmat.Sxx_noiselevel;
%             Sxx_h_0{i_N,i_T} = curmat.Sxx_h_0;
%             Sxx_omega_0{i_N,i_T} = curmat.Sxx_omega_0;
%             Sxx_gamma{i_N,i_T} = curmat.Sxx_gamma;
%             Sxx_fitob{i_N,i_T} = curmat.Sxx_fitob;
%             Syy_noiselevel{i_N,i_T} = curmat.Syy_noiselevel;
%             Syy_h_0{i_N,i_T} = curmat.Syy_h_0;
%             Syy_omega_0{i_N,i_T} = curmat.Syy_omega_0;
%             Syy_gamma{i_N,i_T} = curmat.Syy_gamma;
%             Syy_fitob{i_N,i_T} = curmat.Syy_fitob;
%             Smperpmperp_noiselevel{i_N,i_T} = curmat.Smperpmperp_noiselevel;
%             Smperpmperp_h_0{i_N,i_T} = curmat.Smperpmperp_h_0;
%             Smperpmperp_omega_0{i_N,i_T} = curmat.Smperpmperp_omega_0;
%             Smperpmperp_gamma{i_N,i_T} = curmat.Smperpmperp_gamma;
%             Smperpmperp_fitob{i_N,i_T} = curmat.Smperpmperp_fitob;
%             Smparmpar_noiselevel{i_N,i_T} = curmat.Smparmpar_noiselevel;
%             Smparmpar_h_0{i_N,i_T} = curmat.Smparmpar_h_0;
%             Smparmpar_omega_0{i_N,i_T} = curmat.Smparmpar_omega_0;
%             Smparmpar_gamma{i_N,i_T} = curmat.Smparmpar_gamma;
%             Smparmpar_fitob{i_N,i_T} = curmat.Smparmpar_fitob;
%             Stt_noiselevel{i_N,i_T} = curmat.Stt_noiselevel;
%             Stt_h_0{i_N,i_T} = curmat.Stt_h_0;
%             Stt_omega_0{i_N,i_T} = curmat.Stt_omega_0;
%             Stt_gamma{i_N,i_T} = curmat.Stt_gamma;
%             Stt_fitob{i_N,i_T} = curmat.Stt_fitob;
%             Sww_noiselevel{i_N,i_T} = curmat.Sww_noiselevel;
%             Sww_h_0{i_N,i_T} = curmat.Sww_h_0;
%             Sww_omega_0{i_N,i_T} = curmat.Sww_omega_0;
%             Sww_gamma{i_N,i_T} = curmat.Sww_gamma;
%             Sww_fitob{i_N,i_T} = curmat.Sww_fitob;
% 
%             qbin{i_N,i_T}=curmat.qbin;
%             om_vals{i_N,i_T}=curmat.om_vals;
%             averaging_times{i_N,i_T}=curmat.averaging_times;

        end
    end
save(collectfilepath, "Sxx_noiselevel", "Sxx_h_0", "Sxx_omega_0", "Sxx_gamma", ...
    "Syy_noiselevel", "Syy_h_0", "Syy_omega_0", "Syy_gamma", ...
    "Smperpmperp_noiselevel", "Smperpmperp_h_0", "Smperpmperp_omega_0", "Smperpmperp_gamma", ...
    "Smparmpar_noiselevel", "Smparmpar_h_0", "Smparmpar_omega_0", "Smparmpar_gamma", ...
    "Stt_noiselevel", "Stt_h_0", "Stt_omega_0", "Stt_gamma", ...
    "Sww_noiselevel", "Sww_h_0", "Sww_omega_0", "Sww_gamma");
%     save(collectfilepath,'Sxx', 'Syy', 'Smperpmperp', 'Smparmpar', 'Stt', 'Sww', ...
%         'Sxx_std', 'Syy_std', 'Smperpmperp_std', 'Smparmpar_std', 'Stt_std', 'Sww_std', ...
%         'Sxx_noiselevel', 'Sxx_h_0', 'Sxx_omega_0', 'Sxx_gamma', 'Sxx_fitob', 'Syy_noiselevel', 'Syy_h_0', 'Syy_omega_0', 'Syy_gamma', 'Syy_fitob', ...
%         'Smperpmperp_noiselevel', 'Smperpmperp_h_0', 'Smperpmperp_omega_0', 'Smperpmperp_gamma', 'Smperpmperp_fitob', 'Smparmpar_noiselevel', 'Smparmpar_h_0', 'Smparmpar_omega_0', 'Smparmpar_gamma', 'Smparmpar_fitob', ...
%         'Stt_noiselevel', 'Stt_h_0', 'Stt_omega_0', 'Stt_gamma', 'Stt_fitob', 'Sww_noiselevel', 'Sww_h_0', 'Sww_omega_0', 'Sww_gamma', 'Sww_fitob', ...
%         'om_vals','qbin','runmax','averaging_times');
end