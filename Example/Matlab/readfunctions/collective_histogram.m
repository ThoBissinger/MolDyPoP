%% 0 Initialization
rootdir = "/data/scc/thobi/211201_LongerTime/mxy_3.00";
sampfilename='sampling_output_Dynamics';
collectfilename="samp_Dynamics_Mhistogram";
T_dirs = [".11" ".14" ".165" ".17" ".175" ".18" ".185"];
sqrtN_vals = [16 32 64 128];
runmax_vals = [1000 1000 500 500];

%% 1 For loop
N_N = numel(sqrtN_vals);
N_T = numel(T_dirs);
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    n_run = runmax_vals(i_N);
    for i_T = 1:N_T
        T_str = T_dirs(i_T);
        basedir = sprintf('%s/sqrtN_%d/T_%s',rootdir,sqrtN,T_str);
        M_histogram(basedir,sampfilename,collectfilename,n_run)
    end
end