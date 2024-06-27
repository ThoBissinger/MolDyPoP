function collect_helicity(pathbase,runmax,sim_id,system,L,N,kT)
%COLLECT_GR Collects the two-point correlation fct over all runs in a path.
%   Variables
%   pathbase      Basic directory containing all the run folders, e.g.
%                 /data/scc/thobi/200513_MXYModel/rho_3.00/sqrtN_256/T_.19
%                 /data/scc/thobi/200213_XYModel/lf0/sqrtN_256/T_1.58
%   runmax        maximal number of runs in a directory. 
%   sim_id        the id of the run, used to define
%                 snapshot_SIM_ID_final.out etc.
%   system        kind of system: only mxy or vm makes true sense.
    addpath("/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/ScanningSnapshots")
%%  Initialization of variables
    collectfilepath=sprintf('%s/samp_%s_helicity',pathbase,sim_id);
    collectfilepath_mat=sprintf('%s.mat',collectfilepath);
    snapfilename = sprintf('snapshot_%s_final',sim_id);
    sampfilename = sprintf('sampling_output_%s_helicity',sim_id);
    
    old_size = 0;       % Number of runs already included in collectfile
    if ( isfile(collectfilepath_mat))
        S=load(collectfilepath_mat,'runmax');
        old_size = S.runmax;
    end
%     old_size = 0; % tremove later
    if (old_size < runmax)
%         'H_x','H_y','I_x','I_y','H_s','H_r'
        H_x_collect = zeros(runmax,1);
        H_y_collect = zeros(runmax,1);
        I_x_collect = zeros(runmax,1);
        I_y_collect = zeros(runmax,1);
        H_s_collect = zeros(runmax,1);
        H_r_collect = zeros(runmax,1);
        
        %% Running the loop
        %  Extracts data if no .mat-file is present. Otherwise, collects the
        %  data.
        for runnr = 1:runmax
            curpath=sprintf('%s/run_%i/output',pathbase,runnr);
            matfilepath=sprintf('%s/%s.mat',curpath,sampfilename);
            snapfilepath=sprintf('%s/%s.out',curpath,snapfilename);
            
            disp(curpath);
            if (~ isfile(matfilepath) )
                calculate_helicity(snapfilepath, matfilepath, system, L);
            end
    
            S=load(matfilepath);
            H_x_collect(runnr) = S.H_x;
            H_y_collect(runnr) = S.H_y;
            I_x_collect(runnr) = S.I_x;
            I_y_collect(runnr) = S.I_y;
            H_s_collect(runnr) = S.H_s;
            H_r_collect(runnr) = S.H_r;
            
        end
        H_x = mean(H_x_collect);
        H_y = mean(H_y_collect);
        I_x = mean(I_x_collect);
        I_y = mean(I_y_collect);
        H_r = mean(H_r_collect);
        H_s = mean(H_s_collect);
        I_x_2 = mean(I_x_collect.^2);
        I_y_2 = mean(I_y_collect.^2);

        H_x_std = std(H_x_collect);
        H_y_std = std(H_y_collect);
        I_x_std = std(I_x_collect);
        I_y_std = std(I_y_collect);
        H_r_std = std(H_r_collect);
        H_s_std = std(H_s_collect);
        I_x_2_std = std(I_x_collect.^2);
        I_y_2_std = std(I_y_collect.^2);

        rho_s = 1/2/N*(-(H_x + H_y) - 1/kT * (I_x_2 + I_y_2));


        disp(collectfilepath_mat);
        save(collectfilepath_mat, 'H_x_collect', 'H_y_collect', ...
            'I_x_collect', 'I_y_collect', 'H_r_collect', 'H_s_collect', ...
            'H_x', 'H_y', 'I_x', 'I_y', 'H_r', 'H_s', ...
            'I_x_2', 'I_y_2',...
            'H_x_std', 'H_y_std', 'I_x_std', 'I_y_std', 'H_r_std', 'H_s_std',...
            'I_x_2_std','I_y_2_std',...
            'rho_s','N','kT','runmax');
    
%         disp(collectfilepath_reduced);
%         save(collectfilepath_reduced,'rbin','gr');
    end
end