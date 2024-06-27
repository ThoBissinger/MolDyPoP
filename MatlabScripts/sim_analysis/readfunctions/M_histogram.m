function M_histogram(basedir,sampfilename,collectfilename,n_run)
%M_HISTOGRAM Summary of this function goes here
%   Detailed explanation goes here
    collectfile=sprintf("%s/%s.mat",basedir,collectfilename);
    if ~isfile(collectfile)
        edges=-8:.1:3.4;    
        hist_collect = zeros(1,numel(edges)-1);
        mean_collect = 0;
        var_collect = 0;
        for i_run=1:n_run
            samppath=sprintf('%s/run_%d/output/%s.mat',basedir,i_run,sampfilename);
            fprintf("M_histogram    %s\n",samppath);
            M_cur=load(samppath,"absM").absM;
            [counts, histedges] = histcounts((M_cur-mean(M_cur))/sqrt(var(M_cur)),edges);
            hist_collect = hist_collect + counts;
            mean_collect = mean_collect + mean(M_cur);
            var_collect = var_collect + var(M_cur);
        end
        M_histcounts = hist_collect ./ n_run;
        M_edges=edges(1:end-1);
        M_pdf=M_histcounts/sum(M_histcounts)/(M_edges(2)-M_edges(1));
        M_mean = mean_collect / n_run;
        M_var = var_collect / n_run;
    
        collectfile=sprintf("%s/%s",basedir,collectfilename);
        disp(collectfile);
        save(collectfile,"M_edges","M_histcounts","M_pdf", "M_mean", "M_var");
    else
        fprintf("collectfile already exists, %s\n",collectfile)
    end
end

