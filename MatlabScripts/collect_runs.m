function collect_runs(pathbase,set_specifier,model) %,runmax)
%COLLECT_RUNS Collects a given set of variables over all runs in a path.
%   Variables:
%   pathbase: Basic directory containing all the run folders, e.g.
%       /data/scc/thobi/200513_MXYModel/rho_3.00/sqrtN_256/T_.19
%       /data/scc/thobi/200213_XYModel/lf0/sqrtN_256/T_1.58
%   set_specifier: Sampling set to be analyzed, e.g.
%       'reduced','qreduced',...
%   model: model type (determines which variables are checked). Options:
%       'xy', 'mxy', 'vm' (TODO)
%%%   runmax (redundant): maximal number of runs in a directory. Can be set to 0 or less
%%%       to be auto-determined.

%%  Initialization   
    names=name_vector(set_specifier,model);
%   mfilename:    name of file to be run (should be a .m file)
%   matfilename:  name of matlab file to be generated from runfile
%   collectfile:  file that the collected data will be stored in
    if ( set_specifier == 'eq' || set_specifier == 'integ' )
        mfilename = sprintf('sampling_output_%s.m',set_specifier);
    elseif ( set_specifier == 'eq_samp' )
%         set_specifier = 'eq';
        mfilename = sprintf('eq_sampling_output_eq.m',set_specifier);
    else
        mfilename = sprintf('sampling_output_samp_%s.m',set_specifier);
    end
    matfilename = sprintf('%sat',mfilename);
    collectfile = sprintf('%s/samp_%s.mat',pathbase,set_specifier);

%   Counts the number of subfolders in the pathdir. Should be equal to the 
%   number of runs.
    pathdir = dir(pathbase); 
    runmax=sum([pathdir(~ismember({pathdir.name},{'.','..'})).isdir]);

    mfiles=convertCharsToStrings(sprintf('%s/run_',pathbase)) + (1:runmax) + convertCharsToStrings(sprintf('/output/%s',mfilename));
    matfiles=convertCharsToStrings(sprintf('%s/run_',pathbase)) + (1:runmax) + convertCharsToStrings(sprintf('/output/%s',matfilename));
    if ( set_specifier == 'eq_samp' )
        mfile_aux=convertCharsToStrings(sprintf('%s/run_',pathbase)) + (1:runmax) + convertCharsToStrings(sprintf('/output/eq_sampling_output_eq_1.m'));
        mfileruns=isfile(mfile_aux); % Has ones for all runs that have an mfile (logical array)
    else
        mfileruns=isfile(mfiles); % Has ones for all runs that have an mfile (logical array)
    end
    matfileruns=isfile(matfiles); % Has ones for all runs that already have a matfile (logical array)
    

    total_runnr = sum(mfileruns); % Total number of finished runs.
    
    runmin=find(matfileruns,1) + 1;
    if(isempty(runmin))
        runmin=1;
    end
%     if (runmax <= 0)
%         curfile=sprintf('%s/run_%i/output/%s',pathbase,runmax+1,runfilename);
%         disp(curfile)
%         while(isfile(curfile))
%             runmax=runmax+1;
%             curfile=sprintf('%s/run_%i/output/%s',pathbase,runmax+1,runfilename);
%         end
%     end    
    if (runmax > 0)
%         runmin=1;
%         if (isfile(collectfile))
%             load(collectfile,"averaging_times_collect");
%             runmin = length(averaging_times_collect(:,1)) + 1;
%             clear H_av_collect;
%             fprintf('runmin = %d, runmax = %d\n', runmin,runmax);
%         end
        curpath=sprintf('%s/run_%i/output',pathbase,1);
        disp(curpath);
        
        matfile=sprintf('%s/%s',curpath,matfilename);
        if (~matfileruns(1))
            if ( set_specifier == 'eq_samp' )
                read_eq_output(curpath);
            else
                extract_variables_to_mat(curpath,mfilename,matfilename);
            end
        end
        for name_ind = 1:length(names)
            if (runmin == 1 || ~isfile(collectfile))
%                 curvar = sprintf('%s_cur',names(name_ind));
                curvar = sprintf('%s',names(name_ind));
                load(matfile,curvar);
                
                evalstr=sprintf('%s_collect = zeros(runmax,length(%s));',names(name_ind),curvar);
                eval(evalstr);
                evalstr=sprintf('%s_collect(1,:) = %s;',names(name_ind),curvar);
                eval(evalstr);
            else
                load(collectfile,sprintf('%s_collect',names(name_ind)));
%                 evalstr=sprintf('%s_collect_new = %s_collect;',names(name_ind),names(name_ind));
%                 eval(sprintf('%s_collect_new = %s_collect;',names(name_ind),names(name_ind)));
%                 eval(sprintf('%s_collect = zeros(runmax,length(%s_collect_new(1,:)));',names(name_ind),names(name_ind)));
%                 eval(sprintf('%s_collect(find(matfileruns),:) = %s_collect_new;',names(name_ind),names(name_ind)));
%                 eval(sprintf('clear %s_collect_new;',names(name_ind)));
                eval(sprintf('%s_collect = [%s_collect;zeros(runmax-length(%s_collect(:,1)),length(%s_collect(1,:)))];',names(name_ind),names(name_ind),names(name_ind),names(name_ind)));
            end
        end

    
        for runnr = find(averaging_times_collect(:,end) == 0)'
            curpath=sprintf('%s/run_%i/output',pathbase,runnr);
            matfile=sprintf('%s/%s',curpath,matfilename);
            if ((mfileruns(runnr)))
                disp(curpath);
                 if (~matfileruns(runnr))
                     if ( set_specifier == 'eq_samp' )
                        read_eq_output(curpath);
                     else
                        extract_variables_to_mat(curpath,mfilename,matfilename);
                     end
                 end
                for name_ind = 1:length(names)
%                     curvar = sprintf('%s_cur',names(name_ind));
                    curvar = sprintf('%s',names(name_ind));
                    load(matfile,curvar);
                    evalstr=sprintf('%s_collect(%i,:) = %s;',names(name_ind),runnr,curvar);
    %                 disp(evalstr);
                    eval(evalstr);
                end
            else
                fprintf('No file in %s\n', curpath);
            end
                
        end
    end
%     while(isfile(curpath))
%         disp(curpath);
%         evalstr=sprintf('assignin(''base'', ''%s_cur'',%s)',names(i),names(i));
%         eval(evalstr);
%         runmax=runmax+1;
%         curpath=sprintf('%s/run_%i/output/%s',pathbase,runmax,filename);
%         extract_variables(pathname,filename,names);
%     end
%     finalfilename=sprintf('runcollect_%s',filename);
    disp(collectfile);
    for name_ind = 1:length(names)
        evalstr=sprintf('%s = sum(%s_collect) / total_runnr;',names(name_ind),names(name_ind));
%         evalstr=sprintf('%s = mean(%s_collect(find(%s_collect ~= 0)));',names(name_ind),names(name_ind),names(name_ind));
        eval(evalstr);
%         evalstr=sprintf('clear %s_cur;',names(name_ind));
%         eval(evalstr);
%         save(collectfilename,variables,'-append')
    end
    
    
%     collectfile=sprintf('%s/%s',pathbase,collectfile);
    clear runnr runmax pathbase names name_ind matfilename matfile runnr ...
         curfile curpath curvar evalstr filename ...
         matfileruns mfileruns mfiles matfiles;
    
    save(collectfile);
end
