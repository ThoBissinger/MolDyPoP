function extract_variables_to_mat(pathname,filename,matfilename)
%EXTRACT VARIABLES TO MAT Stores 'names' variables in pathname/filename
%   to their corrseponding cur variables in the .mat-file 
%   pathname/matfilenamein the corresponding folder.
    file=fullfile(pathname,filename);
    run(file);
    matfile=fullfile(pathname,matfilename);
    if (contains(pathname,"/fmxy"))
        H_av = mean(H);
        H_2_av = mean(H_2);
        Hkin_2_av = mean(Hkin_2);
        Hint_2_av = mean(Hint_2);
        H_var = H_2_av - H_av^2;
        temperature_av = mean(temperature);
        temperature_squared_av = mean(temperature_squared);
        temperature_var = temperature_squared_av - temperature_av^2;
        C_0 = H_var / temperature_av^2;
        W_av = mean(W);
        W_2_av = mean(W_2);
        W_var = W_2_av - W_av^2;
    end
    save(matfile);
end
