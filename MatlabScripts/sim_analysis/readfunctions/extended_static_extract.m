function extended_static_extract(pathbase,runmax,runsep,eq_Tmax,sampfile,collectfile,varnames)
%EXTENDED_STATIC_EXTRACT Collects static variables from all run under path,
%with extended time (since simulation is picked up)
%   Variables
%   pathbase      Basic directory containing all the run folders, e.g.
%                 /data/scc/thobi/200513_MXYModel/rho_3.00/sqrtN_256/T_.19
%                 /data/scc/thobi/200213_XYModel/lf0/sqrtN_256/T_1.58
%   runmax        maximal number of runs in a directory. 
%   runmax        number of runs after which a simulation is continued in
%                 the next run (125, in our case)
%   eq_Tmax       equilibration time (must be added to runs)
%   sampfile      name of the sampfile, e.g. sampling_output_Dynamics (no 
%                 .m or .mat suffix). In each pathbase/run/output/... 
%                 folder
%   collectfile   the file into which all the data is to be collected in 
%                 the directory pathbase
%   varnames      name of the variables

%%  Initialization of variables
    collectfilepath=sprintf('%s/%s',pathbase,collectfile);
    collectfilepath_mat=sprintf('%s/%s.mat',pathbase,collectfile);
    old_size = 0;       % Number of runs already included in collectfile
    if ( isfile(collectfilepath_mat))
        load(collectfilepath,'averaging_times_collect');
        old_size = numel(averaging_times_collect(:,1));
    end
    if (old_size < runmax)
        curpath=sprintf('%s/run_1/output',pathbase);
        mfilename=sprintf('%s.m',sampfile);
        matfilename=sprintf('%s.mat',sampfile);
        matfilepath=sprintf('%s/%s.mat',curpath,sampfile);
        if (~ isfile(matfilepath) )
            fprintf('ERROR in extended_static_extract\n      no .mat file under %s\n      abort calculation\n',matfilename);
            return
        end

    %     extract_variables_to_mat(curpath,mfilename,matfilename);
        
        if (old_size == 0)
            load(matfilepath,'M');
            Mx_collect = [M(1,:); zeros(runmax-1,numel(M(1,:)))];
            My_collect = [M(2,:); zeros(runmax-1,numel(M(2,:)))];
            for name_ind = 1:length(varnames)
                load(matfilepath,varnames(name_ind));
                eval(sprintf('%s_collect = [%s;zeros(runmax-1,length(%s))];',varnames(name_ind),varnames(name_ind),varnames(name_ind)));
            end
            old_size = 1;
        else
            for name_ind = 1:length(varnames)
                load(collectfilepath,sprintf('%s_collect',varnames(name_ind)));
                eval(sprintf('%s_collect = [%s_collect;zeros(runmax-old_size,length(%s_collect(1,:)))];',varnames(name_ind),varnames(name_ind),varnames(name_ind)));
            end
        end
    %% Running the loop
    %  Extracts data if no .mat-file is present. Otherwise, collects the
    %  data.
        for runnr = old_size+1:runmax
            curpath=sprintf('%s/run_%i/output',pathbase,runnr);
            mfilepath=sprintf('%s/%s.m',curpath,sampfile);
            matfilepath=sprintf('%s/%s.mat',curpath,sampfile);

            
            if (~ isfile(matfilepath) )
                fprintf('ERROR in extended_static_extract\n      no .mat file under %s\n      abort calculation\n',matfilepath);
                return
            end

            disp(curpath);
            for name_ind = 1:length(varnames)
                load(matfilepath,varnames(name_ind));
                evalstr=sprintf('%s_collect(%i,:) = %s;',varnames(name_ind),runnr,varnames(name_ind));
                eval(evalstr);
            end

        end

        disp(collectfilepath);
        for name_ind = 1:length(varnames)
            evalstr=sprintf('%s = sum(%s_collect) / runmax;',varnames(name_ind),varnames(name_ind));
            eval(evalstr);
        end


        clear runnr runmax pathbase names name_ind matfilename matfile runnr ...
             old_size curfile curpath curvar evalstr filename ...
             matfileruns matfilepath mfileruns mfiles matfiles;
    
        save(collectfilepath);
    end
end
