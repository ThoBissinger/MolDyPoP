%% 0 Initializing data
clear all
close all

sim_id = "eq";
runmax=500;
system = "xy_s";
addpath("/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/StructureFactor")
addpath("/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/ScanningSnapshots")
T_str = [".10" ".20" ".30" ".40" ".50" ".60" ".70" ".80" ".85" ".87" ".89" ".90" ".91" ".93" ".95" ".97" "1.00" "1.03" "1.06" "1.09" "1.10" "1.12" "1.15" "1.18" "1.20" "1.21" "1.24" "1.30" "1.40" "1.50" "1.60" "1.70" "1.80" "1.90" "2.00"];
T_vals = str2double(T_str);
sqrtN_vals = [16,32,64,128,256];
sqrtN_vals = [128,256];
N_N =numel(sqrtN_vals);
N_T =numel(T_vals);
jfunc=@(r) 1;
ufunc=@(r) 0;
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    N = sqrtN^2;
    L = sqrtN;
    for i_T = 1:N_T
        kT = T_vals(i_T);
        if sqrtN == 16
            pathbase = sprintf("/data/scc/thobi/201207_equilibration/xy_s/anneal/sqrtN_%d/T_%s",sqrtN,T_str{i_T});
        else
            pathbase = sprintf("/data/scc/thobi/201207_equilibration/xy_s/scale/sqrtN_%d/T_%s",sqrtN,T_str{i_T});
        end
        fprintf("%s\n",pathbase);
        collect_helicity(pathbase,runmax,sim_id,system,L,N,kT);
    end
end