function collect_chi(pathbase,runmax,sim_id,system,L,qmax)
%COLLECT_CHI Collects the susceptibility fct over all runs in a path (final configuration).
%   Variables
%   pathbase      Basic directory containing all the run folders, e.g.
%                 /data/scc/thobi/200513_MXYModel/rho_3.00/sqrtN_256/T_.19
%                 /data/scc/thobi/200213_XYModel/lf0/sqrtN_256/T_1.58
%   runmax        maximal number of runs in a directory. 
%   sim_id        the id of the run, used to define
%                 snapshot_SIM_ID_final.out etc.
%   system        kind of system. "mxy", "xy", "xy_s", "vm", "fmxy"
%   L             system size
%   qmax          maximum value for q

%%  Initialization of variables
    collectfilepath_reduced=sprintf('%s/samp_%s_chi',pathbase,sim_id);
    collectfilepath_full=sprintf('%s/samp_%s_chi_collect',pathbase,sim_id);
    collectfilepath_full_mat=sprintf('%s.mat',collectfilepath_full);
    snapfilename = sprintf('snapshot_%s_final',sim_id);
    sampfilename = sprintf('sampling_output_%s_chi',sim_id);
    
    old_size = 0;       % Number of runs already included in collectfile
    if ( isfile(collectfilepath_full_mat))
        load(collectfilepath_full_mat,'chi_collect');
        old_size = numel(chi_collect(:,1));
    end
    if (old_size < runmax)
        if (old_size == 0)
            curpath=sprintf('%s/run_1/output',pathbase);
            matfilepath=sprintf('%s/%s.mat',curpath,sampfilename);
            snapfilepath=sprintf('%s/%s.out',curpath,snapfilename);
            disp(curpath);
            if (~ isfile(matfilepath) )
                calculate_chi(snapfilepath, matfilepath, system,L,qmax);
            end
            
            load(matfilepath,'q_vals','chi','chimpar','chimperp', ...
                'chi_mfree','chimpar_mfree','StructFac');
            n_q=numel(q_vals);
            chi_collect = zeros(runmax,n_q);
            chimpar_collect = zeros(runmax,n_q);
            chimperp_collect = zeros(runmax,n_q);
            chi_mfree_collect = zeros(runmax,n_q);
            chimpar_mfree_collect = zeros(runmax,n_q);
            StructFac_collect = zeros(runmax,n_q);
        else
            load(collectfilepath_full_mat,'q_vals',...
                'chi_collect','chimpar_collect','chimperp_collect', ...
                'chi_mfree_collect','chimpar_mfree_collect','StructFac_collect');
            n_q=numel(q_vals);
            chi_collect = [chi_collect;zeros(runmax-old_size,n_q)];
            chimpar_collect = [chimpar_collect;zeros(runmax-old_size,n_q)];
            chimperp_collect = [chimperp_collect;zeros(runmax-old_size,n_q)];
            chi_mfree_collect = [chi_mfree_collect;zeros(runmax-old_size,n_q)];
            chimpar_mfree_collect = [chimpar_mfree_collect;zeros(runmax-old_size,n_q)];
            StructFac_collect = [StructFac_collect;zeros(runmax-old_size,n_q)];
        end
    else
        load(collectfilepath_full_mat,'q_vals',...
            'chi_collect','chimpar_collect','chimperp_collect', ...
            'chi_mfree_collect','chimpar_mfree_collect','StructFac_collect');
        n_q=numel(q_vals);
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
            calculate_chi(snapfilepath, matfilepath, system,L,qmax);
        end

        
        load(matfilepath,'q_vals','chi','chimpar','chimperp', ...
            'chi_mfree','chimpar_mfree','StructFac');
        chi_collect(runnr,:) = chi;
        chimpar_collect(runnr,:) = chimpar;
        chimperp_collect(runnr,:) = chimperp;
        chi_mfree_collect(runnr,:) = chi_mfree;
        chimpar_mfree_collect(runnr,:) = chimpar_mfree;
        StructFac_collect(runnr,:) = StructFac;
    end

    chi=mean(chi_collect);
    chi_err=std(chi_collect)/sqrt(runmax);
    chimpar=mean(chimpar_collect);
    chimpar_err=std(chimpar_collect)/sqrt(runmax);
    chimperp=mean(chimperp_collect);
    chimperp_err=std(chimperp_collect)/sqrt(runmax);
    chi_mfree=mean(chi_mfree_collect);
    chi_mfree_err=std(chi_mfree_collect)/sqrt(runmax);
    chimpar_mfree=mean(chimpar_mfree_collect);
    chimpar_mfree_err=std(chimpar_mfree_collect)/sqrt(runmax);
    StructFac=mean(StructFac_collect);
    StructFac_err=std(StructFac_collect)/sqrt(runmax);

    disp(collectfilepath_full);
    save(collectfilepath_full_mat,'q_vals',...
        'chi_collect','chimpar_collect','chimperp_collect', ...
        'chi_mfree_collect','chimpar_mfree_collect','StructFac_collect');
    
    disp(collectfilepath_reduced);
    save(collectfilepath_reduced,'q_vals','chi','chimpar','chimperp', ...
        'chi_mfree','chimpar_mfree','StructFac', ...
        'chi_err','chimpar_err','chimperp_err', ...
        'chi_mfree_err','chimpar_mfree_err','StructFac_err');
end
