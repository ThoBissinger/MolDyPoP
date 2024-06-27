addpath '/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/StructureFactor'
addpath '/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/ScanningSnapshots'

clear
close all
curpath="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/BondOrder";
plotfold=sprintf("%s/honey_tria_test",curpath);
T_str=[".01" ".09" ".14" "triangle" "honeycomb"];
N_T=numel(T_str)-2;
T_vals=str2double(T_str(1:N_T));
sqrtN=128;
L=74;
rho=sqrtN^2/L^2;
a_tria = sqrt(2/sqrt(3)/(128/74)^2);
a_honey = sqrt(2)*a_tria;
K=24;
% L=74 * [1;sqrt(3)/2];

% rho=sqrtN^2/L(1)/L(2);
pathbase='/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00';
snapname='snapshot_Dynamics_final.out';
grfilename='samp_Dynamics_grSCF';

columnwidth_cm=8.5;

sumfunc=@(r,rbin,g) 2*pi*rho*(rbin(2)-rbin(1))*sum(rbin(rbin<r).*g(rbin<r));

%% Data collection
% cursnapfile=sprintf('%s/sqrtN_%d/T_%s/run_1/output/%s',pathbase,sqrtN,T_str(1),snapname);
% [r_tria,~,~,~] = mxy_snapshot_extract(cursnapfile,'r','xy','t');
% r_honey=[r_tria,r_tria+[0;1/sqrt(3)]];
% r_tria=r_tria*L(1)/sqrtN;
% r_honey=r_honey*sqrt(2)*L(1)/sqrtN;
% r_honey=r_honey(:,r_honey(1,:)<=L(1));
% r_honey=r_honey(:,r_honey(2,:)<=L(2));

[r_tria,L_tria] = trigonal_lattice(L,a_tria);
[r_aux,L_honey] = trigonal_lattice(L,a_honey);
r_honey=[r_aux,r_aux+[0;1/sqrt(3)*a_honey]];
rho_tria = numel(r_tria)/2/L_tria(1)/L_tria(2);
rho_honey= numel(r_honey)/2/L_honey(1)/L_honey(2);

%% Neighbor distance histogram
[r_hon_new,~] = extend_pbc(r_honey,[],L_honey,5);
[Ind_NN,D_NN]=knnsearch(r_hon_new',r_honey','K',K,'Distance','Euclidean');
rmax=1.1*max(max(D_NN));
dr=.01;
r_histbin=0:dr:rmax;
r_hist=zeros(K,numel(r_histbin)-1);
for i_K = 1:K
    r_hist(i_K,:)=histcounts( D_NN(:,i_K),r_histbin') / size(D_NN,1);
end
r_histbin=r_histbin(1:end-1);

%% Plot of distance histogram

c_map=turbo(K);
figure
colormap(turbo(K))
for i_K = 1:K
    plot(r_histbin,r_hist(i_K,:),...
        'Color',c_map(i_K,:));
    hold on;
end
c=colorbar("Ticks",0:K-1);
hx=xlabel('$r_{ij}$','interpreter','latex');
hy=ylabel('$P(r_{ij})$','interpreter','latex');
caxis([0 K-1]);
c.Label.String = '$n$th nearest neighbor';
c.Label.Interpreter = 'latex';
% c.Limits = [0 K-1];
%% 
rmax_honey=[.5 .6 .7,.8,.9,1.00];
rmax_tria=[.3 .4 .5 .6 .7 .8];
n_r_h=numel(rmax_honey);
n_r_t=numel(rmax_tria);
rand_vals=[0 .1 .15 .2 .25];
n_rand=numel(rand_vals);
% rmax_Psi_vals=[-1,0,1.2];
dPsi=.01;
Psi_histbin=-dPsi/2:dPsi:(1+dPsi/2);
n_hist=numel(Psi_histbin);
N_histbin=0:25;
k=6;
Psi_hist_hon=zeros(n_rand,n_r_h,n_hist-1);
N_hist_hon=zeros(n_rand,n_r_h,numel(N_histbin)-1);

Psi_hist_tria=zeros(n_rand,n_r_t,n_hist-1);
N_hist_tria=zeros(n_rand,n_r_t,numel(N_histbin)-1);
for i_rand=1:n_rand
    for i_r=1:n_r_h
        [psi_k,N_N] = bondorder(r_honey + rand_vals(i_rand)*rand(size(r_honey)),L_honey,6,rmax_honey(i_r));
        aux_hist=histcounts( abs(psi_k),Psi_histbin');
        Psi_hist_hon(i_rand,i_r,:)=reshape(aux_hist/(numel(r_honey)/2)/(Psi_histbin(2)-Psi_histbin(1)),[],1);
        aux_hist=histcounts( N_N,N_histbin');
        N_hist_hon(i_rand,i_r,:)=reshape(aux_hist,[],1)/(numel(r_honey)/2);
    end
    for i_r=1:n_r_t
        [psi_k,N_N] = bondorder(r_tria + rand_vals(i_rand)*rand(size(r_tria)),L_tria,6,rmax_honey(i_r));
        aux_hist=histcounts( abs(psi_k),Psi_histbin');
        Psi_hist_tria(i_rand,i_r,:)=reshape(aux_hist/(numel(r_tria)/2)/(Psi_histbin(2)-Psi_histbin(1)),[],1);
        aux_hist=histcounts( N_N,N_histbin');
        N_hist_tria(i_rand,i_r,:)=reshape(aux_hist,[],1)/(numel(r_tria)/2);
    end
end
Psi_histbin=.5*(Psi_histbin(2:end)+Psi_histbin(1:end-1));
N_histbin=N_histbin(1:end-1);

%%
c_map=lines(100);
% c_map=c_map([35 55 70 85 95],:);
for i_rand=1:n_rand
    figure
    for i_r=1:n_r_h
        dispname=sprintf("$r_{\\textrm{max}}=%.2f$",rmax_honey(i_r));
        plot(Psi_histbin,reshape(Psi_hist_hon(i_rand,i_r,:),size(Psi_histbin)),'-*',...
            'DisplayName',dispname,...
            'Color',c_map(i_r,:),...
            'LineWidth',1.2);
        hold on;
    end
    xlabel("$|\psi_6|$",'interpreter','latex');
    ylabel("$P(|\psi_6|)$",'interpreter','latex');
    hleg{i_r}=legend('NumColumns',2,'Location','northeast','interpreter','latex');
    set(gca,'FontName','cmr14');
    hleg{i_r}=legend('NumColumns',2,'Location','northeast','interpreter','latex');
    title(sprintf("honeycomb, $N = %d$, rand $= %.2f$",sqrtN,rand_vals(i_rand)),"interpreter",'latex');
%     ylim([0 14]);
% 
%     set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
%     figname=sprintf("%s/honey_hist_randfix_r_%.2f",plotfold,rmax_honey(i_r));
%     fprintf("Printing %s.png\n",figname);
%     exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
end

%%
c_map=turbo(100);
c_map=c_map([35 55 70 85 95],:);
for i_rand=1:n_rand
    figure
    for i_r=1:n_r_h
        dispname=sprintf("$r_{\\textrm{max}}=%.2f$",rmax_honey(i_r));
        plot(Psi_histbin,reshape(Psi_hist_tria(i_rand,i_r,:),size(Psi_histbin)),'-*',...
            'DisplayName',dispname,...
            'Color',c_map(i_r,:),...
            'LineWidth',1.2);
        hold on;
    end
    xlabel("$|\psi_6|$",'interpreter','latex');
    ylabel("$P(|\psi_6|)$",'interpreter','latex');
    hleg{i_r}=legend('NumColumns',2,'Location','northeast','interpreter','latex');
    set(gca,'FontName','cmr14');
    hleg{i_r}=legend('NumColumns',2,'Location','northeast','interpreter','latex');
    title(sprintf("triangular, $N = %d$, rand $= %.2f$",sqrtN,rand_vals(i_rand)),"interpreter",'latex');
    ylim([0 14]);

    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    figname=sprintf("%s/tria_hist_randfix_r_%.2f",plotfold,rmax_honey(i_r));
    fprintf("Printing %s.png\n",figname);
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
end


%%
c_map=turbo(100);
c_map=c_map([35 55 70 85 95],:);

for i_r=1:n_r_h
    figure
    for i_rand=1:n_rand
        dispname=sprintf("rand $= %.2f$",rand_vals(i_rand));
        plot(Psi_histbin,reshape(Psi_hist_hon(i_rand,i_r,:),size(Psi_histbin)),'-*',...
            'DisplayName',dispname,...
            'Color',c_map(i_rand,:),...
            'LineWidth',1.2);
        hold on;
    end
    xlabel("$|\psi_6|$",'interpreter','latex');
    ylabel("$P(|\psi_6|)$",'interpreter','latex');
    hleg{i_r}=legend('NumColumns',2,'Location','northeast','interpreter','latex');
    set(gca,'FontName','cmr14');
    hleg{i_r}=legend('NumColumns',2,'Location','northeast','interpreter','latex');
    title(sprintf("honeycomb, $N = %d, r_{\\textrm{max}}= %.2f$",sqrtN,rmax_honey(i_r)),"interpreter",'latex');
    ylim([0 14]);

    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    figname=sprintf("%s/honey_hist_rfix_rand_%.2f",plotfold,rand_vals(i_rand));
    fprintf("Printing %s.png\n",figname);
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
end

%%
c_map=turbo(100);
c_map=c_map([35 55 70 85 95],:);

for i_r=1:n_r_h
    figure
    for i_rand=1:n_rand
        dispname=sprintf("rand $= %.2f$",rand_vals(i_rand));
        plot(Psi_histbin,reshape(Psi_hist_tria(i_rand,i_r,:),size(Psi_histbin)),'-*',...
            'DisplayName',dispname,...
            'Color',c_map(i_rand,:),...
            'LineWidth',1.2);
        hold on;
    end
    xlabel("$|\psi_6|$",'interpreter','latex');
    ylabel("$P(|\psi_6|)$",'interpreter','latex');
    hleg{i_r}=legend('NumColumns',2,'Location','northeast','interpreter','latex');
    set(gca,'FontName','cmr14');
    title(sprintf("triangular, $N = %d, r_{\\textrm{max}}= %.2f$",sqrtN,rmax_honey(i_r)),"interpreter",'latex');
    ylim([0 14]);

    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    figname=sprintf("%s/tria_hist_rfix_rand_%.2f",plotfold,rand_vals(i_rand));
    fprintf("Printing %s.png\n",figname);
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
end