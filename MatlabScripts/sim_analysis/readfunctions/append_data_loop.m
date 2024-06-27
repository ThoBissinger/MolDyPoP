%% 0 Initializing data
clear all
close all
%  Includes paths to other data.
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/FitFuncs'));

%  section_select chooses the task to be performed.
%  section_select = 1
%       Loops through all folders with data and converts .m files into .mat
%       files. Collects data for all runs of a sqrtN-T-pair at the level of 
%       the sqrtN/T/... folder
%  section_select = 2
%       Combines the collected data into a single file with a lot of
%       information

section_select = 1:2;
% section_select = 2;

%  Some parameters that are relevant for the functions.
%  dataset_id       The dataset ID switch gives different modes for this 
%                   program to run. Can be used to have different data sets
%                   for analysis etc.
%                   Values:
%                   eq_mxy          For equilibration analysis (MXY model)
%                   eq_xy           For equilibration analysis (XY model)
%                   dynamics_mxy    For dynamics analysis  (MXY)
%                   dynamics_mxy_fullT  As above, but with more
%                                   temperatures
%                   dynamics_mxy_fullT  As above, but with better
%                                   resolution of q variables
%                   LepriRuffo_mxy  For Lepri&Ruffo low temperature spin
%                                   wave dynamics simulation (MXY)
%                   LepriRuffo_xy   For Lepri&Ruffo low temperature spin
%                                   wave dynamics simulation (XY)
%                   
%  model            The simulation model. 'mxy', 'xy', 'vm', 'fvm
%  rho              Density. Relevant for mxy and vm
%  sampfilename     Name of the sampling file to be read
%  collectfilename  Name of the file that the collected data should be 
%                   stored in for each temperature (so the file at the
%                   level of sqrtN/T/...) 
%  extractfile      The file that the final collection of data should be
%                   collected into
%  dirs             Root directory / directories.
%  sqrtN_dirs       Subdirectory name for sqrtN values
%  T_dirs_cell      Cell array with T subdirectory names. If it is left
%                   empty, then all subdirectories of the sqrtN_dirs will
%                   be used
%  runmax           Number of runs. Required for data extraction.
%  varnames         Variable names: all variables that should be gathered.
%                   Easily generated with the name_vector function.
%  calc_crit_temp_switch     Checks whether the critical temperatures
%                   should be calculated from the magnetization during
%                   after data extraction.
% dataset_id='LepriRuffo_mxy';
% dataset_id='LepriRuffo_fmxy';
% dataset_id='LepriRuffo_xy';
% dataset_id='LepriRuffo_extended_mxy';
% dataset_id='LepriRuffo_extended_fmxy';
dataset_id='dynamics_fmxy_AdjustedTime';
% dataset_id='dynamics_mxy_LinearTime';
% dataset_id='dynamics_xy_LinearTime';
% dataset_id='dynamics_fmxy_LinearTime';
% dataset_id='static_mxy';
% dataset_id='dynamics_mxy';
% dataset_id='dynamics_mxy_fullT';
% dataset_id='dynamics_mxy_better_q';
% dataset_id='eq_mxy';
% dataset_id='eq_xy';
% dataset_id='eq_xy_s';

if (strcmp(dataset_id,'LepriRuffo_mxy'))
    model='mxy';
    rho=3.00;
    sampfilename='sampling_output_LepriRuffo';
    collectfilename='samp_LepriRuffo';
    extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_LepriRuffo.mat";
    dirs=["/data/scc/thobi/210727_LepriRuffo_GoodResolution/mxy_3.00"];
%     dirs=["/data/scc/thobi/210528_LepriRuffoRerun/mxy_3.00"];
    sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128", "sqrtN_256"];
%     sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128"];
    L_vals=[9.25,18.5,37,74,148];
    T_dirs_cell={};
    T_dirs_cell={"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20"};
    runmax=250;
    varnames=name_vector(dataset_id);
    calc_crit_temp_switch=0;
elseif (strcmp(dataset_id,'LepriRuffo_fmxy'))
    model='fmxy';
    sampfilename='sampling_output_LepriRuffo';
    collectfilename='samp_LepriRuffo';
    extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/fmxy/fmxy_LepriRuffo.mat";
    dirs=["/data/scc/thobi/210727_LepriRuffo_GoodResolution/fmxy"];
%     dirs=["/data/scc/thobi/210528_LepriRuffoRerun/mxy_3.00"];
    sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128", "sqrtN_256"];
%     sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128"];
    L_vals=[9.25,18.5,37,74,148];
    T_dirs_cell={};
%     T_dirs_cell={"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.185" "T_.19" "T_.195" "T_.20" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
    runmax=250;
    varnames=name_vector(dataset_id);
    calc_crit_temp_switch=0;
elseif (strcmp(dataset_id,'LepriRuffo_xy'))
    model='xy';
    rho=3.00;
    sampfilename='sampling_output_LepriRuffo';
    collectfilename='samp_LepriRuffo';
    extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/xy/xy_LepriRuffo.mat";
%     dirs=["/data/scc/thobi/210528_LepriRuffoRerun/xy"];
    dirs=["/data/scc/thobi/210727_LepriRuffo_GoodResolution/xy"];
    sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128", "sqrtN_256"];
    L_vals=[16,32,64,128,256];
%     T_dirs_cell={"T_.10", "T_.20", "T_.30", "T_.40", "T_.50", "T_.60", "T_.70", "T_.80", "T_.90", "T_1.00"};
    T_dirs_cell={};
    runmax=125;
    varnames=name_vector(dataset_id);
    calc_crit_temp_switch=0;
elseif (strcmp(dataset_id,'LepriRuffo_extended_mxy') || strcmp(dataset_id,'LepriRuffo_extended_fmxy'))
    if ( strcmp(dataset_id,'LepriRuffo_extended_mxy') ) 
        model='mxy';
        rho=3.00;
        extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_LepriRuffo_extended.mat";
        dirs=["/data/scc/thobi/211126_ShortTime_GoodResolution/mxy_3.00"];
    elseif ( strcmp(dataset_id,'LepriRuffo_extended_fmxy') )
        extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/fmxy/fmxy_LepriRuffo_extended.mat";
        dirs=["/data/scc/thobi/211126_ShortTime_GoodResolution/fmxy"];
    end
    sampfilename='sampling_output_Dynamics';
    collectfilename='samp_Dynamics';
    
%     dirs=["/data/scc/thobi/210528_LepriRuffoRerun/mxy_3.00"];
    sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128", "sqrtN_256"];
%     sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128"];
    L_vals=[9.25,18.5,37,74,148];
    T_dirs_cell={};
%     T_dirs_cell={"T_.07" "T_.09" "T_.11" "T_.14" "T_.15" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.20"};
    runmax=250;
    varnames=name_vector(dataset_id);
    calc_crit_temp_switch=0;
elseif (strcmp(dataset_id,'eq_mxy'))
    model='mxy';
    rho=3.00;
    sampfilename='sampling_output_eq';
    collectfilename='samp_eq';
    extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_eq.mat";
    dirs=["/data/scc/thobi/201207_equilibration/mxy_3.00"];
    sqrtN_dirs=["anneal/sqrtN_16" "scale/sqrtN_32", "scale/sqrtN_64", "scale/sqrtN_128" "scale/sqrtN_256"];
%     sqrtN_dirs=["scale/sqrtN_256"];
    L_vals=[9.25,18.5,37,74,148];
    T_dirs_cell={};

    runmax=125;
    varnames=name_vector(dataset_id);
    calc_crit_temp_switch=1;
    
elseif (strcmp(dataset_id,'eq_xy') || strcmp(dataset_id,'eq_xy_s'))
    model='xy';
    rho=3.00;
    sampfilename='sampling_output_eq';
    collectfilename='samp_eq';
    if ( strcmp(dataset_id,'eq_xy') )
        extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/xy/eq_xy.mat";
        dirs=["/data/scc/thobi/201207_equilibration/xy"];
        sqrtN_dirs=["anneal/sqrtN_16" "anneal/sqrtN_32", "anneal/sqrtN_64", "anneal/sqrtN_128" "scale/sqrtN_256"];
        L_vals=[16,32,64,128,256];
    else
        extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/xy/eq_xy.mat";
        dirs=["/data/scc/thobi/201207_equilibration/xy_s"];
%         sqrtN_dirs=["anneal/sqrtN_16" "scale/sqrtN_32", "scale/sqrtN_64", "scale/sqrtN_128" "scale/sqrtN_256"];
%         L_vals=[16,32,64,128,256];
        sqrtN_dirs=["anneal/sqrtN_16" "scale/sqrtN_32", "scale/sqrtN_64", "scale/sqrtN_128"];
        L_vals=[16,32,64,128];
    end
    T_dirs_cell={};
    runmax=125;
    varnames=name_vector('eq_xy');
    calc_crit_temp_switch=1;
    
elseif (strcmp(dataset_id,'dynamics_mxy') || strcmp(dataset_id,'dynamics_mxy_fullT') || strcmp(dataset_id,'dynamics_mxy_better_q'))
    model='mxy';
    rho=3.00;
    sampfilename='sampling_output_Dynamics';
    collectfilename='samp_integ';
    extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_dynamics.mat";  
    dirs=["/data/scc/thobi/210606_DynamicsReRerun/mxy_3.00"];
%     dirs=["/data/scc/thobi/210317_Dynamics/mxy_3.00"];
    sqrtN_dirs=["sqrtN_256"];
    L_vals=[148];
    T_dirs_cell={"T_.11", "T_.17", "T_.18", "T_.185", "T_.19", "T_.195", "T_.20", "T_.24"};
%     T_dirs_cell={};
    runmax=375;
%     varnames=name_vector(dataset_id);
    varnames=name_vector("dynamics_mxy");
    calc_crit_temp_switch=0;
    if (strcmp(dataset_id,'dynamics_mxy_fullT'))
        runmax=125;
%         T_dirs_cell={"T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
        T_dirs_cell={"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
        collectfilename='samp_integ_fullT';
        extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_dynamics_fullT.mat";  
        calc_crit_temp_switch=1;
    elseif (strcmp(dataset_id,'dynamics_mxy_better_q'))
        dirs=["/data/scc/thobi/210703_DynamicsBetterQbin/mxy_3.00"];
        runmax=122;
%         T_dirs_cell={"T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
        T_dirs_cell={};
        collectfilename='samp_dynamics';
        extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_dynamics_better_q.mat";  
        calc_crit_temp_switch=0;
    end
    
elseif ( strcmp(dataset_id,'dynamics_mxy_LinearTime') || strcmp(dataset_id,'static_mxy') )
    model='mxy';
    rho=3.00;
    sampfilename='sampling_output_Dynamics';
    if strcmp(dataset_id,'dynamics_mxy_LinearTime')
        collectfilename='samp_Dynamics';
        extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_dynamics_LinearTime.mat";  
    elseif strcmp(dataset_id,'static_mxy')
        collectfilename='samp_static';
        extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_static_LinearTime.mat";  
    end
    dirs=["/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00"];
    sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128", "sqrtN_256"];
    L_vals=[9.25,18.5,37,74,148];
%     sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128"];
    T_dirs_cell={};
%     sqrtN_dirs=["sqrtN_128", "sqrtN_256"];
%     T_dirs_cell={"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
    runmax=500;
    varnames=name_vector(dataset_id);
    calc_crit_temp_switch=1;
elseif (strcmp(dataset_id,'dynamics_fmxy_LinearTime') || strcmp(dataset_id,'static_fmxy'))
    model='fmxy';
    rho=3.00;
    sampfilename='sampling_output_Dynamics';
    if strcmp(dataset_id,'dynamics_fmxy_LinearTime')
        collectfilename='samp_Dynamics';
        extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/fmxy/fmxy_dynamics_LinearTime.mat"; 
    elseif strcmp(dataset_id,'static_fmxy')
        collectfilename='samp_static';
        extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/fmxy/fmxy_static_LinearTime.mat"; 
    end
    dirs=["/data/scc/thobi/210715_LinearTimeSampling/fmxy"];
    sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128", "sqrtN_256"];
    L_vals=[9.25,18.5,37,74,148];
%     sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128"];
    T_dirs_cell={};
%     sqrtN_dirs=["sqrtN_128", "sqrtN_256"];
%     T_dirs_cell={"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
    runmax=500;
    varnames=name_vector(dataset_id);
    calc_crit_temp_switch=1;
elseif (strcmp(dataset_id,'dynamics_fmxy_AdjustedTime') )
    model='fmxy';
    rho=3.00;
    sampfilename='sampling_output_Dynamics';
    collectfilename='samp_Dynamics';
    extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/fmxy/fmxy_dynamics_AdjustedTime.mat"; 
    dirs=["/data/scc/thobi/211201_LongerTime/fmxy"];
    sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128"];
    L_vals=[9.25,18.5,37,74 ];
    T_dirs_cell={};
    T_dirs_cell={"T_.11" "T_.14" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185"};
%     sqrtN_dirs=["sqrtN_128", "sqrtN_256"];
%     T_dirs_cell={"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
    runmax_vals=[3000,1500,375,125,0];
    varnames=name_vector(dataset_id);
    calc_crit_temp_switch=1;
elseif (strcmp(dataset_id,'dynamics_xy_LinearTime') || strcmp(dataset_id,'static_xy'))
    model='xy';
    rho=3.00;
    sampfilename='sampling_output_Dynamics';
    if strcmp(dataset_id,'dynamics_xy_LinearTime')
        collectfilename='samp_Dynamics';
        extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/xy/xy_dynamics_LinearTime.mat";  
    elseif strcmp(dataset_id,'static_xy')
        collectfilename='samp_static';
        extractfile="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/xy/xy_static_LinearTime.mat"; 
    end
    dirs=["/data/scc/thobi/210715_LinearTimeSampling/xy"];
    sqrtN_dirs=["sqrtN_16" "sqrtN_32", "sqrtN_64", "sqrtN_128", "sqrtN_256"];
%     sqrtN_dirs=["sqrtN_256"];
%     T_dirs_cell={"T_1.66"};
    L_vals=[16,32,64,128,256];
    T_dirs_cell={};
    runmax=250;
    varnames=name_vector(dataset_id);
    calc_crit_temp_switch=1;
end


%% 1 Reading all files if necessary
%  Creates collectfiles for each parameter combination. That is, in this
%  part the .m-files in sqrtN/T/runnr/output/... are converted into
%  .mat-files (if necessary) and the variables specified in varnames are
%  extracted from the .mat-files and collected into
%  sqrtN/T/collectfilename. There is also averaging being performed. This
%  is mainly just the application of collect_runs in each sqrtN/T
%  directory.
section_ind = 1;
if (ismember(section_ind,section_select))
    for i=1:numel(dirs)
        fprintf('%s\n',dirs(i))
        for i_N = 3:numel(sqrtN_dirs)
            runmax=runmax_vals(i_N);
            diffdir = dir(fullfile(dirs(i),sqrtN_dirs(i_N),'*')); % improve by specifying the file extension.
            if (isempty(T_dirs_cell))
                T_dirs = setdiff({diffdir([diffdir.isdir]).name},{'.','..'});
            else
                T_dirs = T_dirs_cell;
            end
            for i_T = 1:numel(T_dirs)
                for i_run = 1:runmax
                    curpath=sprintf('%s/%s/%s/run_%d',dirs(i),sqrtN_dirs(i_N),T_dirs{i_T},i_run);
                    curmat=sprintf('%s/output/%s.mat',curpath,sampfilename);
                    disp(curmat);
                    load(curmat,"absM");
                    ACF_q0_absM=absM(1)*absM;
                    save(curmat,"ACF_q0_absM",'-append');
                end
            end
        end
    end
end
