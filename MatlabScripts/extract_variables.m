function extract_variables(pathname,filename,names)
%EXTRACT VARIABLES Stores variables in pathname/filename
%   to their corrseponding cur variables in the base workspace.
%   The chosen variables are given to the function in the vector
%   of strings names
    file=fullfile(pathname,filename);
    run(file);
    for i = 1:length(names)
%         evalstr=['assignin(\'base\'',names(i)+"cur",',',names(i),')']
       evalstr=sprintf('assignin(''base'', ''%s_cur'',%s)',names(i),names(i));
       eval(evalstr);
%         assignin('base',names(i)+"cur",N)
    end
end