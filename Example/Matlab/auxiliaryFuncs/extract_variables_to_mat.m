function extract_variables_to_mat(pathname,filename,matfilename)
%EXTRACT VARIABLES TO MAT Stores 'names' variables in pathname/filename
%   to their corrseponding cur variables in the .mat-file 
%   pathname/matfilenamein the corresponding folder.
    file=fullfile(pathname,filename);
    run(file);
    matfile=fullfile(pathname,matfilename);
%    for i = 1:length(names)
%       evalstr=sprintf('%s_cur = %s;',names(i),names(i));
%       eval(evalstr);
%       evalstr=sprintf('clear %s;',names(i));
%       eval(evalstr);
%    end
    save(matfile);
end
