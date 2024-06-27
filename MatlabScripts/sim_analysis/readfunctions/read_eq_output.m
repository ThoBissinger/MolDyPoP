function read_eq_output(eqfold)
%READ_EQ_OUTPUT Reads output after sampled equilibration
%   Gathers and reads everything in a given folder
    addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
    addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
    if not(isfolder(eqfold))
        fprintf('no folder %s\n',eqfold)
        return
    end
    oldfold=cd(eqfold);
    final_mfile='eq_sampling_output_eq.m';
%     if (~ isfile(sprintf('%s/%s',eqfold,final_mfile)))
    if (~ isfile(final_mfile) )
        i=0;
        while (isfile(sprintf('eq_sampling_output_eq_%d.m',i+1)))
            i=i+1;
        end
        if (i ~= 0)
            mfile=sprintf('eq_sampling_output_eq_%d.m',i);
            copyfile(mfile,final_mfile);
        else
            fprintf('ERROR: no eq sampling output in %s\n',eqfold);
            return
        end
    end
%     if (~ isfile(sprintf('%s/%s',eqfold,mfile)))
    matfile=sprintf('%sat',final_mfile);
    if (~ isfile(matfile) )
        extract_variables_to_mat(eqfold,final_mfile,matfile)
        load(matfile)
        if (~exist('averaging_times','var'))
            fprintf('ERROR: Incorrect output. Reason could be a segfault. Directory %s\n', eqfold)
            return
        end
        zero_entries=find(averaging_times == 0);
        zero_entries(end+1) = numel(averaging_times) + 1;
        for k = 2:numel(zero_entries)-1
            averaging_times(zero_entries(k):zero_entries(k+1) - 1) = ...
                averaging_times(zero_entries(k):zero_entries(k+1) - 1) + ...
                2 * averaging_times(zero_entries(k)-1) - averaging_times(zero_entries(k)-2);
            H(k-1) = [];
            absM(k-1) = [];
            H_2(k-1) = [];
            Hint_2(k-1) = [];
            Hkin_2(k-1) = [];
            M(:,k-1) = [];
            M_2(k-1) = [];
            M_4(k-1) = [];
            temperature(k-1) = [];
            temperature_squared(k-1) = [];
            Theta(k-1) = [];
            W(k-1) = [];
            W_2(k-1) = [];
        end
        M_x = M(1,:);
        M_y = M(2,:);
        save(matfile)
%         for i = floor(numel(averaging_times))/20:numel(averaging_times)
%             
%         end
    end
%     delete(matfile)
    cd(oldfold)
end


