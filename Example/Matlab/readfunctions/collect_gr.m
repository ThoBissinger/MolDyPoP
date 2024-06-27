function collect_gr(pathbase,runmax,sim_id,system)
%COLLECT_GR Collects the two-point correlation fct over all runs in a path.
%   Variables
%   pathbase      Basic directory containing all the run folders, e.g.
%                 /data/scc/thobi/200513_MXYModel/rho_3.00/sqrtN_256/T_.19
%                 /data/scc/thobi/200213_XYModel/lf0/sqrtN_256/T_1.58
%   runmax        maximal number of runs in a directory. 
%   sim_id        the id of the run, used to define
%                 snapshot_SIM_ID_final.out etc.
%   system        kind of system: only mxy or vm makes true sense.

%%  Initialization of variables
    collectfilepath_reduced=sprintf('%s/samp_%s_gr',pathbase,sim_id);
    collectfilepath_full=sprintf('%s/samp_%s_gr_collect',pathbase,sim_id);
    collectfilepath_full_mat=sprintf('%s.mat',collectfilepath_full);
    snapfilename = sprintf('snapshot_%s_final',sim_id);
    sampfilename = sprintf('sampling_output_%s_gr',sim_id);
    
    old_size = 0;       % Number of runs already included in collectfile
    if ( isfile(collectfilepath_full_mat))
        load(collectfilepath_full_mat,'gr_collect','rbin_minlength');
        old_size = numel(gr_collect(:,1));
    end
    if (old_size < runmax)
%         curpath=sprintf('%s/run_1/output',pathbase);
%         matfilepath=sprintf('%s/%s.mat',curpath,sampfilename);
%         snapfilepath=sprintf('%s/%s.out',curpath,snapfilename);
%         if (~ isfile(matfilepath) )
%             calculate_gr(snapfilepath, matfilepath, system);
%         end

    %     extract_variables_to_mat(curpath,mfilename,matfilename);
        
        if (old_size == 0)
            gr_collect = cell(runmax,1);
            rbin_collect = cell(runmax,1);
            rbin_minlength = Inf;
        else
            load(collectfilepath_full_mat,'gr_collect','rbin_collect','rbin_minlength');
            gr_collect = {gr_collect;cell(runmax-old_size,1)};
            rbin_collect = {rbin_collect;cell(runmax-old_size,1)};
        end
    else
        load(collectfilepath_full_mat,'gr_collect','rbin_collect','rbin_minlength');
    end
    %% Running the loop
    %  Extracts data if no .mat-file is present. Otherwise, collects the
    %  data.
    for runnr = old_size+1:runmax
        curpath=sprintf('%s/run_%i/output',pathbase,runnr);
        matfilepath=sprintf('%s/%s.mat',curpath,sampfilename);
        snapfilepath=sprintf('%s/%s.out',curpath,snapfilename);
        
        disp(curpath);
        if (~ isfile(matfilepath) )
            calculate_gr(snapfilepath, matfilepath, system);
        end

        
        load(matfilepath,'gr','rbin');
        gr_collect{runnr} = gr;
        rbin_collect{runnr} = rbin;
%         rbin_minlength = min(numel(rbin),rbin_minlength);
    end

    rbin_minlength=Inf;
    for runnr = 1:runmax
        rbin_minlength=min(numel(rbin_collect{runnr}),rbin_minlength);
    end
    gr_vec=zeros(runmax,rbin_minlength);
    rbin = rbin_collect{1}(1:rbin_minlength);
    
    for runnr = 1:runmax
        gr_vec(runnr,:) = gr_collect{runnr}(1:rbin_minlength);
    end
    gr = mean(gr_vec);

    disp(collectfilepath_full);
    save(collectfilepath_full,'rbin','gr',...
        "rbin_collect","gr_collect","gr_vec");

    disp(collectfilepath_reduced);
    save(collectfilepath_reduced,'rbin','gr');
end
