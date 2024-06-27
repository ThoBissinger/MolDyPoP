function M_TimeEvolution(basedir,sampfilename,collectfilename,n_run)
%M_HISTOGRAM Summary of this function goes here
%   Detailed explanation goes here
    
    collectfile=sprintf("%s/%s.mat",basedir,collectfilename);
%     if ~isfile(collectfile) 
    if 1 == 1
        curfile = sprintf("%s/run_1/output/%s",basedir,sampfilename);
        load(curfile,'averaging_times');
        n_t = numel(averaging_times);
        ACF_q0_M_collect = zeros(n_run,n_t);
        angles = zeros(n_run,n_t);
        Mx_collect = zeros(n_run,n_t);
        My_collect = zeros(n_run,n_t);
        for i_run = 1:n_run
            fprintf("%-4d  %s\n", i_run, basedir);
            curfile = sprintf("%s/run_%d/output/%s",basedir,i_run,sampfilename);
            load(curfile,'M');
            Mx_collect(i_run,:) = M(1,:);
            My_collect(i_run,:) = M(2,:);
            ACF_q0_M_collect(i_run,:) = M(1,1) * M(1,:) + M(2,1) * M(2,:);
            angle_init = atan2(M(2,1),M(1,1));
            Mx = cos(-angle_init) * M(1,:) - sin(-angle_init) * M(2,:);
            My = cos(-angle_init) * M(2,:) + sin(-angle_init) * M(1,:);
            angles(i_run,:) = atan2(My,Mx);
        end
        phi_diffs=[zeros(n_run,1),angles(:,2:end)-angles(:,1:end-1)];
        phi_diffs(phi_diffs>pi)=phi_diffs(phi_diffs>pi) - 2*pi;
        phi_diffs(phi_diffs<-pi)=phi_diffs(phi_diffs<-pi) + 2*pi;
        phi_rel_to_0 = cumsum(phi_diffs,2);
        absM_collect = sqrt(Mx_collect.^2 + My_collect.^2);
        absM_av = mean(mean(absM_collect));
        Mx_av = mean(Mx_collect);
        My_av = mean(My_collect);
        ACF_q0_M = mean(ACF_q0_M_collect);
        M_ang_MSD = mean(phi_rel_to_0.^2);
    
        disp(collectfile);
        save(collectfile,...
            "averaging_times","angles","Mx_collect","My_collect",...
            "Mx_av","My_av","M_ang_MSD",...
            "ACF_q0_M", "ACF_q0_M_collect",...
            "absM_collect", "absM_av",...
            "phi_diffs","phi_rel_to_0");
    else
        fprintf("collectfile already exists, %s\n",collectfile)
    end
end

