addpath '/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/StructureFactor'
addpath '/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/ScanningSnapshots'

clear
close all
curpath="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/BondOrder";
plotfold=sprintf("%s/honey_tria_test",curpath);
T_str=[".01" ".09" ".14"];
N_T=numel(T_str);
T_vals=str2double(T_str(1:N_T));
sqrtN=128;
L=74;
rho=sqrtN^2/L^2;
K=24;
% L=74 * [1;sqrt(3)/2];

% rho=sqrtN^2/L(1)/L(2);
pathbase='/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00';
snapname='snapshot_Dynamics_final.out';
grfilename='samp_Dynamics_grSCF';

columnwidth_cm=8.5;

sumfunc=@(r,rbin,g) 2*pi*rho*(rbin(2)-rbin(1))*sum(rbin(rbin<r).*g(rbin<r));

%% Data collection
for i_T = 1:N_T
    curgrfile=sprintf('%s/sqrtN_%d/T_%s/%s.mat',pathbase,sqrtN,T_str(i_T),grfilename);
    cursnapfile=sprintf('%s/sqrtN_%d/T_%s/run_4/output/%s',pathbase,sqrtN,T_str(i_T),snapname);
    [r,~,~,~] = mxy_snapshot_extract(cursnapfile,'r','mxy','t');
    S{i_T}.r=r;
end

%% Neighbor distance histogram
for i_T = 1:N_T
    [r_new,~] = extend_pbc(S{i_T}.r,[],L,5);
    [Ind_NN,D_NN]=knnsearch(r_new',S{i_T}.r','K',K,'Distance','Euclidean');
    rmax=1.1*max(max(D_NN));
    dr=.01;
    r_histbin=0:dr:rmax;
    r_hist=zeros(K,numel(r_histbin)-1);
    for i_K = 1:K
        r_hist(i_K,:)=histcounts( D_NN(:,i_K),r_histbin') / size(D_NN,1) / dr;
    end
    r_histbin=r_histbin(1:end-1);
    S{i_T}.r_histbin=r_histbin;
    S{i_T}.r_hist=r_hist;
end
    
%% Plot of distance histogram
for i_T = 1:N_T
    c_map=turbo(K);
    figure
    colormap(turbo(K))
    for i_K = 2:K
        plot(S{i_T}.r_histbin,S{i_T}.r_hist(i_K,:),...
            'Color',c_map(i_K,:));
        hold on;
    end
    c=colorbar("Ticks",0:K-1);
    hx=xlabel('$r_{ij}$','interpreter','latex');
    hy=ylabel('$P(r_{ij})$','interpreter','latex');
    caxis([0 K-1]);
    c.Label.String = '$n$th nearest neighbor';
    c.Label.Interpreter = 'latex';

    xline(.66);
    xline(1.06);
    title(sprintf("$N = %d,T = %.2f$",sqrtN,T_vals(i_T)),"interpreter",'latex');
end
