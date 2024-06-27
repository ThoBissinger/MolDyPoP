clear all
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));


dataset_id='eq_mxy';
% dataset_id='eq_NoSpin';
if (strcmp(dataset_id,'eq_mxy'))
    dirs = ["/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_256/T_.01" ...
        "/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_256/T_.09" ...
        "/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_256/T_.15" ...
        "/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_256/T_.19" ...
        "/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_256/T_.25" ...
        "/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_256/T_.31"];
    T_str = ["T_.01" "T_.09" "T_.15" "T_.19" "T_.25" "T_.31"];
    T_vals = [.01 .09 .15 .19 .25 .31];
    runmax=125;
    storefile='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_grselect.mat';
    grfilename='gr_extract_snapshot_eq_final.mat';
    snapfilename='snapshot_eq_final.out';
elseif (strcmp(dataset_id,'eq_NoSpin'))
%     dirs = ["/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.01" ...
%         "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.09" ...
%         "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.15" ...
%         "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.19" ...
%         "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.25" ...
%         "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.31"];
dirs = ["/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.01" ...
        "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.03" ...
        "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.05" ...
        "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.07" ...
        "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.09" ...
        "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.11" ...
        "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.13" ...
        "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.15" ...
        "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.17" ...
        "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.19" ...
        "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.21" ...
        "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.25" ...
        "/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/T_.31"];
    T_str = ["T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.15" "T_.17" "T_.19" "T_.21" "T_.25" "T_.31"];
    T_vals = [.01 .03 .05 .07 .09 .11 .13 .15 .17 .19 .21 .25 .31];
    runmax=50;
    storefile='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/NoSpin/NoSpin_grselect.mat';
    grfilename='gr_extract_snapshot_Dynamics_final.mat';
    snapfilename='snapshot_Dynamics_final.out';
end

% dirs = ["/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_64/T_.01" ...
%     "/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_64/T_.09" ...
%     "/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_64/T_.15" ...
%     "/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_64/T_.19"];

N_dirs = numel(dirs);
% rbin=[];
for i = 1 : N_dirs
%     gr_av{i}=[];
    for i_run = 1 : runmax
        curfile=sprintf('%s/run_%d/output/%s',dirs(i),i_run,snapfilename);
        curgrfile=sprintf('%s/run_%d/output/%s',dirs(i),i_run,grfilename);
        if (isfile(curfile))
            disp(curfile)
            if (~ isfile(curgrfile))
                [r,~,~,~] = mxy_snapshot_extract(curfile,'r','mxy');
                [curgr,currbin,rw]=twopointcorr(r(1,:),r(2,:),.01,100,0);
                save(curgrfile,'curgr','currbin','rw');
            else
                load(curgrfile);
            end
%             if (isempty(gr_av{i}))
%                 gr_av{i} = curgr;
%             else
%                 gr_av{i} = gr_av{i} + curgr;
%             end
            files{i,i_run} = curfile;
            gr{i,i_run} = curgr;
            rbin{i,i_run} = currbin;
        end
    end
%     gr_av{i} = gr_av{i} / runmax;
end

save(storefile, 'files', 'gr', 'rbin','T_str','T_vals',...
    'dirs','runmax','grfilename','snapfilename','dataset_id', 'N_dirs')