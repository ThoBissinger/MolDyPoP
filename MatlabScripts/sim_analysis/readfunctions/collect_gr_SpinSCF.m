function collect_gr_SpinSCF(pathbase,runmax,sim_id,system,L,dr,n_part)
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
    collectfilepath=sprintf('%s/samp_%s_grSCF',pathbase,sim_id);
    collectfilepath_mat=sprintf('%s.mat',collectfilepath);
    snapfilename = sprintf('snapshot_%s_final',sim_id);
    sampfilename = sprintf('sampling_output_%s_grSCF',sim_id);
    
    old_size = 0;       % Number of runs already included in collectfile
    if ( isfile(collectfilepath_mat))
        S=load(collectfilepath_mat,'runmax');
        old_size = S.runmax;
    end
    old_size = 0;
    if (old_size < runmax)
        gr_collect = cell(runmax,1);
        Cm_collect = cell(runmax,1);
        Cm_par_collect = cell(runmax,1);
        Cm_perp_collect = cell(runmax,1);
        Cm_cross_collect = cell(runmax,1);
        Cw_collect = cell(runmax,1);
        Cvx_collect = cell(runmax,1);
        Cvy_collect = cell(runmax,1);

        absM_collect = zeros(runmax,1);
        s_par_2_collect = zeros(runmax,1);
        s_par_4_collect = zeros(runmax,1);
        s_perp_2_collect = zeros(runmax,1);
        s_perp_4_collect = zeros(runmax,1);
        sum_cos_2te_collect = zeros(runmax,1);
        sum_cos_4te_collect = zeros(runmax,1);
        sum_cos_6te_collect = zeros(runmax,1);
        sum_cos_8te_collect = zeros(runmax,1);
        %% Running the loop
        %  Extracts data if no .mat-file is present. Otherwise, collects the
        %  data.
        for runnr = 1:runmax
            curpath=sprintf('%s/run_%i/output',pathbase,runnr);
            matfilepath=sprintf('%s/%s.mat',curpath,sampfilename);
            snapfilepath=sprintf('%s/%s.out',curpath,snapfilename);
            
            disp(curpath);
            if (~ isfile(matfilepath) )
                calculate_gr_spin(snapfilepath, matfilepath, system, L, dr, n_part);
            end
    
            S=load(matfilepath);
            gr_collect{runnr} = S.gr_hist;
            Cm_collect{runnr} = S.Cm;
            Cm_par_collect{runnr} = S.Cm_par;
            Cm_perp_collect{runnr} = S.Cm_perp;
            Cm_cross_collect{runnr} = S.Cm_cross;
            Cw_collect{runnr} = S.Cw;
            Cvx_collect{runnr} = S.Cvx;
            Cvy_collect{runnr} = S.Cvy;
            
            absM_collect(runnr) = S.mnorm;
            s_par_2_collect(runnr) = S.s_par_2;
            s_par_4_collect(runnr) = S.s_par_4;
            s_perp_2_collect(runnr) = S.s_perp_2;
            s_perp_4_collect(runnr) = S.s_perp_4;
            sum_cos_2te_collect(runnr) = S.sum_cos_2te;
            sum_cos_4te_collect(runnr) = S.sum_cos_4te;
            sum_cos_6te_collect(runnr) = S.sum_cos_6te;
            sum_cos_8te_collect(runnr) = S.sum_cos_8te;
    %         rbin_minlength = min(numel(rbin),rbin_minlength);
        end
        rbin = S.rbin;

        gr = mean(cell2mat(gr_collect));
        gr_std = std(cell2mat(gr_collect));
        Cm = mean(cell2mat(Cm_collect));
        Cm_std = std(cell2mat(Cm_collect));
        Cm_par = mean(cell2mat(Cm_par_collect));
        Cm_par_std = std(cell2mat(Cm_par_collect));
        Cm_perp = mean(cell2mat(Cm_perp_collect));
        Cm_perp_std = std(cell2mat(Cm_perp_collect));
        Cm_cross = mean(cell2mat(Cm_cross_collect));
        Cm_cross_std = std(cell2mat(Cm_cross_collect));
        Cw = mean(cell2mat(Cw_collect));
        Cw_std = std(cell2mat(Cw_collect));
        Cvx = mean(cell2mat(Cvx_collect));
        Cvx_std = std(cell2mat(Cvx_collect));
        Cvy = mean(cell2mat(Cvy_collect));
        Cvy_std = std(cell2mat(Cvy_collect));

        Cm_preavg=mean(cell2mat(Cm_collect)./cell2mat(gr_collect));
        Cm_preavg_std=std(cell2mat(Cm_collect)./cell2mat(gr_collect));
        Cm_par_preavg=mean(cell2mat(Cm_par_collect)./cell2mat(gr_collect));
        Cm_par_preavg_std=std(cell2mat(Cm_par_collect)./cell2mat(gr_collect));
        Cm_perp_preavg=mean(cell2mat(Cm_perp_collect)./cell2mat(gr_collect));
        Cm_perp_preavg_std=std(cell2mat(Cm_perp_collect)./cell2mat(gr_collect));
        Cm_cross_preavg=mean(cell2mat(Cm_cross_collect)./cell2mat(gr_collect));
        Cm_cross_preavg_std=std(cell2mat(Cm_cross_collect)./cell2mat(gr_collect));
        Cw_preavg=mean(cell2mat(Cw_collect)./cell2mat(gr_collect));
        Cw_preavg_std=std(cell2mat(Cw_collect)./cell2mat(gr_collect));
        Cvx_preavg=mean(cell2mat(Cvx_collect)./cell2mat(gr_collect));
        Cvx_preavg_std=std(cell2mat(Cvx_collect)./cell2mat(gr_collect));
        Cvy_preavg=mean(cell2mat(Cvy_collect)./cell2mat(gr_collect));
        Cvy_preavg_std=std(cell2mat(Cvy_collect)./cell2mat(gr_collect));

        Cm_preavg(isnan(Cm_preavg))=0;
        Cm_preavg_std(isnan(Cm_preavg_std))=0;
        Cm_par_preavg(isnan(Cm_par_preavg))=0;
        Cm_par_preavg_std(isnan(Cm_par_preavg_std))=0;
        Cm_perp_preavg(isnan(Cm_perp_preavg))=0;
        Cm_perp_preavg_std(isnan(Cm_perp_preavg_std))=0;
        Cm_cross_preavg(isnan(Cm_cross_preavg))=0;
        Cm_cross_preavg_std(isnan(Cm_cross_preavg_std))=0;
        Cw_preavg(isnan(Cw_preavg))=0;
        Cw_preavg_std(isnan(Cw_preavg_std))=0;
        Cvx_preavg(isnan(Cvx_preavg))=0;
        Cvx_preavg_std(isnan(Cvx_preavg_std))=0;
        Cvy_preavg(isnan(Cvy_preavg))=0;
        Cvy_preavg_std(isnan(Cvy_preavg_std))=0;

        gr_absM_2 = mean(absM_collect.^2.*cell2mat(gr_collect));
        gr_absM_2_std = std(absM_collect.^2.*cell2mat(gr_collect));
        
        absM = mean(absM_collect);
        absM_std = std(absM_collect);
        M_2 = mean(absM_collect.^2);
        M_2_std = std(absM_collect.^2);
        M_4 = mean(absM_collect.^4);
        M_4_std = std(absM_collect.^4);

        s_par_2 = mean(s_par_2_collect);
        s_par_2_std = std(s_par_2_collect);
        s_par_4 = mean(s_par_4_collect);
        s_par_4_std = std(s_par_4_collect);
        s_perp_2 = mean(s_perp_2_collect);
        s_perp_2_std = std(s_perp_2_collect);
        s_perp_4 = mean(s_perp_4_collect);
        s_perp_4_std = std(s_perp_4_collect);
        sum_cos_2te = mean(sum_cos_2te_collect);
        sum_cos_2te_std = std(sum_cos_2te_collect);
        sum_cos_4te = mean(sum_cos_4te_collect);
        sum_cos_4te_std = std(sum_cos_4te_collect);
        sum_cos_6te = mean(sum_cos_6te_collect);
        sum_cos_6te_std = std(sum_cos_6te_collect);
        sum_cos_8te = mean(sum_cos_8te_collect);
        sum_cos_8te_std = std(sum_cos_8te_collect);
    
        disp(collectfilepath_mat);
        save(collectfilepath_mat, 'gr', 'Cm', 'Cm_par', 'Cm_perp', ...
            'Cm_cross', 'Cw', 'Cvx', 'Cvy', ...
            'gr_std', 'Cm_std', 'Cm_par_std', 'Cm_perp_std', ...
            'Cm_cross_std', 'Cw_std', 'Cvx_std', 'Cvy_std', ...
            'gr_absM_2','gr_absM_2_std',...
            'Cm_preavg','Cm_preavg_std','Cm_par_preavg','Cm_par_preavg_std',...
            'Cm_perp_preavg','Cm_perp_preavg_std','Cm_cross_preavg','Cm_cross_preavg_std',...
            'Cw_preavg','Cw_preavg_std','Cvx_preavg','Cvx_preavg_std','Cvy_preavg','Cvy_preavg_std',...
            'rbin',...
            'absM_collect','absM','absM_std','M_2','M_2_std',...
            'M_4','M_4_std','dr',...
            'runmax',...
            's_par_2','s_par_2_std','s_par_4','s_par_4_std',....
            's_perp_2','s_perp_2_std','s_perp_4','s_perp_4_std',....
            'sum_cos_2te','sum_cos_2te_std','sum_cos_4te','sum_cos_4te_std',....
            'sum_cos_6te','sum_cos_6te_std','sum_cos_8te','sum_cos_8te_std');
    
%         disp(collectfilepath_reduced);
%         save(collectfilepath_reduced,'rbin','gr');
    end
end