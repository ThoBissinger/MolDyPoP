clear all
close all
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts
% addpath /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/ExternalCodes
addpath(genpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts'));
addpath(genpath('../..'));
addpath(genpath('C:/Users/tbiss/Thesis/MatlabScripts'));
saveswitch=1;
plotpathbase='/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/Snapshots_xy';

runnr=1;

sqrtN = 256;
T_vals = [".10", ".40", ".60", "1.00", "1.50", "1.60", "1.70", "2.00", "2.50", "3.00"];

for i_T = 1:numel(T_vals)
    snapname = "/data/scc/thobi/201207_equilibration/xy/scale/sqrtN_256/T_" + T_vals(i_T) + "/run_" + runnr + "/output/equilibration_final_eq.out";
    single_snap_fig(snapname,'xy',i_T,'spins');

    xlim([0,200]);
    ylim([0,200])



    figname=sprintf('%s/snap_sqrtN_%d_T_%s',plotpathbase,sqrtN,T_vals(i_T));
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
%         exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
end