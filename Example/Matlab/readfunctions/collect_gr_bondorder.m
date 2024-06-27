function collect_gr_bondorder(pathbase,runmax,sim_id,system,L,dr,dPsi,n_part)
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
    collectfilepath=sprintf('%s/samp_%s_bondorder',pathbase,sim_id);
    collectfilepath_mat=sprintf('%s.mat',collectfilepath);
    snapfilename = sprintf('snapshot_%s_final',sim_id);
    sampfilename = sprintf('sampling_output_%s_bondorder',sim_id);
    
    old_size = 0;       % Number of runs already included in collectfile
    if ( isfile(collectfilepath_mat))
        S=load(collectfilepath_mat,'runmax');
        old_size = S.runmax;
    end
    old_size = 0;
    if (old_size < runmax)
        gr_collect = cell(runmax,1);
        Gcorrs_collect = cell(runmax,1);
        Psimean_collect = cell(runmax,1);
        Psiabsmean_collect = cell(runmax,1);
        absPsimean_collect = cell(runmax,1);
        Psihist_collect = cell(runmax,1);
        r_hist_collect = cell(runmax,1);
        
        %% Running the loop
        %  Extracts data if no .mat-file is present. Otherwise, collects the
        %  data.
        for runnr = 1:runmax
            curpath=sprintf('%s/run_%i/output',pathbase,runnr);
            matfilepath=sprintf('%s/%s.mat',curpath,sampfilename);
            snapfilepath=sprintf('%s/%s.out',curpath,snapfilename);
            
            disp(curpath);
%             if (~ isfile(matfilepath) )
                calculate_gr_bondorder(snapfilepath, matfilepath, system, L, dr, dPsi, n_part);
%             end
    
            S=load(matfilepath);
            gr_collect{runnr} = S.gr_hist;
            Gcorrs_collect{runnr} = S.Gcorrs;
            Psimean_collect{runnr} = S.Psimean;
            Psiabsmean_collect{runnr} = S.Psiabsmean;
            absPsimean_collect{runnr} = S.absPsimean;
            Psihist_collect{runnr} = S.Psi_hist;
            r_hist_collect{runnr} = S.r_hist;
            
        end
        rbin = S.rbin;
        Psibin = .5*(S.Psibin(1:end-1) + S.Psibin(2:end));
        r_hist_bin = S.r_hist_bin;
        k_vals = S.k_vals;
        rmax_vals = S.rmax_vals;
        n_k=numel(k_vals);
        n_rmax=numel(rmax_vals);
        n_r=numel(rbin);
        K = size(S.r_hist,1);
        n_r_hist_bin=numel(r_hist_bin);

        gr = mean(cell2mat(gr_collect));
        gr_std = std(cell2mat(gr_collect));
        
        auxGcorr=permute(reshape(cell2mat(Gcorrs_collect),n_k,[],n_rmax,n_r),[2 1 3 4]);
        Gcorrs_mean = reshape(mean(auxGcorr,1),n_k,n_rmax,n_r);
        Gcorrs_abs_mean = reshape(mean(abs(auxGcorr),1),n_k,n_rmax,n_r);
        Gcorrs_std = reshape(std(auxGcorr,1),n_k,n_rmax,n_r);

        Gcorrs_jkn = reshape(jackknife(@mean,auxGcorr),runmax,n_k,n_rmax,n_r);
        Gcorrs_bias=(runmax-1)*(reshape(mean(Gcorrs_jkn),size(Gcorrs_mean)) - Gcorrs_mean);
        Gcorrs_bcmean=Gcorrs_mean-Gcorrs_bias;
        Gcorrs_pseudovals=runmax*reshape(Gcorrs_mean,1,n_k,n_rmax,n_r) - (runmax-1)*Gcorrs_jkn;
        Gcorrs_jackvar=1/(runmax*(runmax-1))*sum((Gcorrs_pseudovals-reshape(Gcorrs_bcmean,1,n_k,n_rmax,n_r)).^2);
        Gcorrs_ste=reshape(sqrt(Gcorrs_jackvar),size(Gcorrs_mean));

        auxPsihist=permute(reshape(cell2mat(Psihist_collect),n_k,[],n_rmax,numel(Psibin)),[2 1 3 4]);
        Psihist_mean = reshape(mean(auxPsihist,1),n_k,n_rmax,numel(Psibin));
        Psihist_std = reshape(std(auxPsihist,1),n_k,n_rmax,numel(Psibin));

        auxPsimean=permute(reshape(cell2mat(Psimean_collect),n_k,[],n_rmax),[2 1 3]);
        Psimean = reshape(mean(auxPsimean,1),n_k,n_rmax);
        Psimean_std = reshape(std(auxPsimean,1),n_k,n_rmax);

        auxPsiabsmean=permute(reshape(cell2mat(Psiabsmean_collect),n_k,[],n_rmax),[2 1 3]);
        Psiabsmean = reshape(mean(auxPsiabsmean,1),n_k,n_rmax);
        Psiabsmean_std = reshape(std(auxPsiabsmean,1),n_k,n_rmax);

        auxabsPsimean=permute(reshape(cell2mat(absPsimean_collect),n_k,[],n_rmax),[2 1 3]);
        absPsimean = reshape(mean(auxabsPsimean,1),n_k,n_rmax);
        absPsimean_std = reshape(std(auxabsPsimean,1),n_k,n_rmax);
        
        auxrhist=permute(reshape(cell2mat(r_hist_collect),K,[],n_r_hist_bin),[2 1 3]);
        r_hist_mean = reshape(mean(auxrhist,1),K,[]);
        r_hist_std = reshape(std(auxrhist,1),K,[]);
    
        disp(collectfilepath_mat);
        save(collectfilepath_mat, 'rbin', 'k_vals', 'rmax_vals', 'Psibin', ...
            'gr', 'Gcorrs_mean', 'Psimean', 'Psiabsmean', 'absPsimean',...
            'Psihist_mean',...
            'Gcorrs_abs_mean','Gcorrs_ste','Gcorrs_bcmean','Gcorrs_jkn',...
            'dr','runmax',...
            'gr_std', 'Gcorrs_std', 'Psimean_std', 'Psiabsmean_std', 'absPsimean_std',...
            'Psihist_std',...
            'r_hist_bin','r_hist_mean','r_hist_std');
    
%         disp(collectfilepath_reduced);
%         save(collectfilepath_reduced,'rbin','gr');
    end
end