function collect_runs(pathbase,runmax,sampfile,collectfile,varnames)
%COLLECT_RUNS Collects a given set of variables over all runs in a path.
%   Variables
%   pathbase      Basic directory containing all the run folders, e.g.
%                 /data/scc/thobi/200513_MXYModel/rho_3.00/sqrtN_256/T_.19
%                 /data/scc/thobi/200213_XYModel/lf0/sqrtN_256/T_1.58
%   runmax        maximal number of runs in a directory. 
%   sampfile      name of the sampfile, e.g. sampling_output_integ (no .m or
%                 .mat suffix). In each pathbase/run/output/... folder
%   collectfile   the file into which all the data is to be collected in 
%                 the directory pathbase
%   varnames      name of the variables

%%  Initialization of variables
    collectfilepath_reduced=sprintf('%s/%s',pathbase,collectfile);
    collectfilepath_full=sprintf('%s/%s_collect',pathbase,collectfile);
    collectfilepath_full_mat=sprintf('%s.mat',collectfilepath_full);
    % Checks whether anything has to be done in this directory
%     readswitch = 1;
    old_size = 0;       % Number of runs already included in collectfile
    if ( isfile(collectfilepath_full_mat))
        load(collectfilepath_full_mat,'averaging_times_collect');
        old_size = numel(averaging_times_collect(:,1));
%         if (numel(averaging_times_collect(:,1)) == runmax)
%             readswitch = 0;
%         end
    end
    if (old_size < runmax)
        curpath=sprintf('%s/run_1/output',pathbase);
        mfilename=sprintf('%s.m',sampfile);
        matfilename=sprintf('%s.mat',sampfile);
        matfilepath=sprintf('%s/%s.mat',curpath,sampfile);
        if (~ isfile(matfilepath) )
            extract_variables_to_mat(curpath,mfilename,matfilename);
        end

    %     extract_variables_to_mat(curpath,mfilename,matfilename);
        
        if (old_size == 0)
            for name_ind = 1:length(varnames)
                if (varnames{name_ind} == "M" || varnames{name_ind} == "Mx" || varnames{name_ind} == "My")
                    load(matfilepath,"M");
                    Mx_collect = [M(1,:); zeros(runmax-1,numel(M(1,:)))];
                    My_collect = [M(2,:); zeros(runmax-1,numel(M(2,:)))];
                elseif (varnames{name_ind} == "P" || varnames{name_ind} == "Px" || varnames{name_ind} == "Py")
                    load(matfilepath,"P");
                    Px_collect = [P(1,:); zeros(runmax-1,numel(P(1,:)))];
                    Py_collect = [P(2,:); zeros(runmax-1,numel(P(2,:)))];
                elseif (varnames{name_ind} == "ACF_q0_absM")
                    load(matfilepath,"absM");
                    ACF_q0_absM = absM(1) * absM;
                    save(matfilepath,"ACF_q0_absM","-append");
                    ACF_q0_absM_collect = [ACF_q0_absM; zeros(runmax-1,numel(ACF_q0_absM))];
                else
                    load(matfilepath,varnames{name_ind});
                    eval(sprintf('%s_collect = [%s;zeros(runmax-1,length(%s))];',varnames{name_ind},varnames{name_ind},varnames{name_ind}));
                end
            end
            old_size = 1;
        else
            for name_ind = 1:length(varnames)
                load(collectfilepath_full,sprintf('%s_collect',varnames{name_ind}));
                eval(sprintf('%s_collect = [%s_collect;zeros(runmax-old_size,length(%s_collect(1,:)))];',varnames{name_ind},varnames{name_ind},varnames{name_ind}));
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
                extract_variables_to_mat(curpath,mfilename,matfilename);
            end

            disp(curpath);
            for name_ind = 1:length(varnames)
                if (varnames{name_ind} == "M" || varnames{name_ind} == "Mx" || varnames{name_ind} == "My")
                    load(matfilepath,"M");
                    Mx_collect(runnr,:) = M(1,:);
                    My_collect(runnr,:) = M(2,:);
                elseif (varnames{name_ind} == "P" || varnames{name_ind} == "Px" || varnames{name_ind} == "Py")
                    load(matfilepath,"P");
                    Px_collect(runnr,:) = P(1,:);
                    Py_collect(runnr,:) = P(2,:);
                elseif (varnames{name_ind} == "ACF_q0_absM")
                    load(matfilepath,"absM");
                    ACF_q0_absM = absM(1) * absM;
                    save(matfilepath,"ACF_q0_absM","-append");
                    ACF_q0_absM_collect(runnr,:)= ACF_q0_absM;
                else
                    load(matfilepath,varnames{name_ind});
                    evalstr=sprintf('%s_collect(%i,:) = %s;',varnames{name_ind},runnr,varnames{name_ind});
                    eval(evalstr);
                end
            end

        end

        disp(collectfilepath_full);
        for name_ind = 1:length(varnames)
            if (varnames{name_ind} == "M" || varnames{name_ind} == "Mx" || varnames{name_ind} == "My")
                Mx = sum(Mx_collect) /runmax;
                My = sum(My_collect) /runmax;
            elseif (varnames{name_ind} == "P" || varnames{name_ind} == "Px" || varnames{name_ind} == "Py")
                Px = sum(Px_collect) /runmax;
                Py = sum(Py_collect) /runmax;
%             elseif (varnames{name_ind} == "ACF_q0_absM")
%                 ACF_q0_absM
            else
                evalstr=sprintf('%s = sum(%s_collect) / runmax;',varnames{name_ind},varnames{name_ind});
                eval(evalstr);
            end
        end


        clear runnr pathbase names name_ind matfilename matfile runnr ...
             old_size curfile curpath curvar evalstr filename ...
             matfileruns matfilepath mfileruns mfiles matfiles;
    
        save(collectfilepath_full);

        disp(collectfilepath_reduced);
        save(collectfilepath_reduced,varnames{:});
    end
end
