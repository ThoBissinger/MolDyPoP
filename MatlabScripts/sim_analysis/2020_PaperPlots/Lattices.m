clear all
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
saveswitch=1;

model="mxy";
% model="NoSpin";

if ( model == "mxy")
    T_str={'T_.01', 'T_.03', 'T_.05', 'T_.09', 'T_.16', 'T_.20'};
    pathbase="/data/scc/thobi/201207_equilibration/mxy_3.00/scale/sqrtN_256";
    pathfinal="run_115/output/equilibration_final_eq.out";
    figbase="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/mxy_snap";
elseif ( model == "NoSpin")
    T_str={'T_.01', 'T_.03', 'T_.05', 'T_.09', 'T_.16', 'T_.20'};
    pathbase="/data/scc/thobi/210625_NoSpinInteraction/mxy_3.00/sqrtN_256/";
    pathfinal="run_1/output/snapshot_Dynamics_final.out";
    figbase="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/NoSpin_snap";
end

for i_T = 1:numel(T_str)
    snapfilename=sprintf('%s/%s/%s',pathbase,T_str{i_T},pathfinal);
    single_snap_fig(snapfilename,'mxy',i_T,'pos')
    figure(i_T)
    xlim([0 15])
    ylim([0 15])
    figname=sprintf('%s_%s',figbase,T_str{i_T});
    if(saveswitch == 1)
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
%         exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
end

