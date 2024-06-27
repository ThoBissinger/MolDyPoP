addpath '/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/StructureFactor'
addpath '/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/ScanningSnapshots'

clear
close all
T_str=[".01" ".09" ".14" "triangle" "honeycomb"];
N_T=numel(T_str)-2;
T_vals=str2double(T_str(1:N_T));
sqrtN=128;
L=74;
L=74 * [1;sqrt(3)/2];
% rho=sqrtN^2/L^2;
rho=sqrtN^2/L(1)/L(2);
pathbase='/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00';
snapname='snapshot_Dynamics_final.out';
grfilename='samp_Dynamics_grSCF';

sumfunc=@(r,rbin,g) 2*pi*rho*(rbin(2)-rbin(1))*sum(rbin(rbin<r).*g(rbin<r));

%% Data collection
S=cell(size(T_str));
for i_T = 1:N_T
    curgrfile=sprintf('%s/sqrtN_%d/T_%s/%s.mat',pathbase,sqrtN,T_str(i_T),grfilename);
    cursnapfile=sprintf('%s/sqrtN_%d/T_%s/run_1/output/%s',pathbase,sqrtN,T_str(i_T),snapname);
    Scur=load(curgrfile,'rbin','gr');
    S{i_T}.rbin=Scur.rbin;
    S{i_T}.gr=Scur.gr;
%     [r,~,~,~] = mxy_snapshot_extract(cursnapfile,'r','mxy','s');
    [r,~,~,~] = mxy_snapshot_extract(cursnapfile,'r','mxy','t');
    r=r*L(1)/sqrtN;
    S{i_T}.r=r;     
end
[r_tria,~,~,~] = mxy_snapshot_extract(cursnapfile,'r','xy','t');
r_honey=[r_tria,r_tria+[0;1/sqrt(3)]];
r_tria=r_tria*L(1)/sqrtN;
r_honey=r_honey*sqrt(2)*L(1)/sqrtN;
r_honey=r_honey(:,r_honey(1,:)<=L(1));
r_honey=r_honey(:,r_honey(2,:)<=L(2));
S{4}.r=r_tria;
S{5}.r=r_honey;

rho=2.99;
%% n_Area calculation
rmax_vals=0:.01:3;
nArea_vals=zeros(N_T,numel(rmax_vals));
for i_T=1:N_T
    for i_r = 1:numel(rmax_vals)
        nArea_vals(i_T,i_r)=sumfunc(rmax_vals(i_r),S{i_T}.rbin,S{i_T}.gr);
    end
end

%% Rho plot
figure
c_map=turbo(100);
c_map=c_map([10,40,90],:);
markers=["s" "o" "^"];
for i_T=1:N_T
    dispname=sprintf("T=%.2f",T_vals(i_T));
    plot(rmax_vals,nArea_vals(i_T,:),'DisplayName',dispname,...
        'Color',c_map(i_T,:),...
        'LineWidth',2,...
        'Marker',markers(i_T),'MarkerIndices',i_T:20:300);
    hold on;
    xlim([.3 1.5]);
end
legend show;
xlabel("$r$","interpreter","latex");
ylabel("$N_g(r)$","interpreter","latex");



%% 
rmax_Psi_vals=[0,.6:.05:1.3];
rmax_Psi_vals=[-1,0,.7,1.05];
% rmax_Psi_vals=[-1,0,1.2];
dPsi=.03;
Psi_histbin=-dPsi/2:dPsi:(1+dPsi/2);
N_histbin=0:25;
k=6;
Psi_hist=zeros(N_T+2,numel(rmax_Psi_vals),numel(Psi_histbin)-1);
N_hist=zeros(N_T+2,numel(rmax_Psi_vals),numel(N_histbin)-1);
for i_T = 1:N_T+2
    for i_r=1:numel(rmax_Psi_vals)
        if rmax_Psi_vals(i_r) < 0
            [psi_k,N_N] = bondorder(S{i_T}.r,L,k,"Delaunay");
        else
            [psi_k,N_N] = bondorder(S{i_T}.r,L,k,rmax_Psi_vals(i_r));
        end
        aux_hist=histcounts( abs(psi_k),Psi_histbin');
        Psi_hist(i_T,i_r,:)=reshape(aux_hist/(numel(r)/2)/(Psi_histbin(2)-Psi_histbin(1)),[],1);
        aux_hist=histcounts( N_N,N_histbin');
        N_hist(i_T,i_r,:)=reshape(aux_hist,[],1)/(numel(r)/2);

    end
end
Psi_histbin=.5*(Psi_histbin(2:end)+Psi_histbin(1:end-1));
N_histbin=N_histbin(1:end-1);

%%
c_map=hsv(numel(rmax_Psi_vals));
for i_T = 1:N_T+2
    figure
    for i_r=3:numel(rmax_Psi_vals)
        dispname=sprintf("$r_{\\textrm{max}}=%.2f$",rmax_Psi_vals(i_r));
        plot(Psi_histbin,reshape(Psi_hist(i_T,i_r,:),size(Psi_histbin)),...
            'DisplayName',dispname,...
            'Color',c_map(i_r,:),...
            'LineWidth',1.2);
        hold on;
    end
    dispname="Delaunay";
    plot(Psi_histbin,reshape(Psi_hist(i_T,1,:),size(Psi_histbin)),...
        'DisplayName',dispname,...
        'Color','k','Marker','^',...
        'LineWidth',1.2);
    dispname="Topological";
    plot(Psi_histbin,reshape(Psi_hist(i_T,2,:),size(Psi_histbin)),...
        'DisplayName',dispname,...
        'Color','k','Marker','s',...
        'LineWidth',1.2);
    
    hleg{i_r}=legend('NumColumns',2,'Location','northeast','interpreter','latex');
    if i_T <= N_T
        title(sprintf("$N = %d,T = %.2f$",sqrtN,T_vals(i_T)),"interpreter",'latex');
    else
        title(sprintf("$N = %d$, %s",sqrtN,T_str(i_T)),"interpreter",'latex');
    end
end


%%
c_map=hsv(numel(rmax_Psi_vals));
for i_T = 1:N_T+2
    
    barmat=zeros(numel(rmax_Psi_vals),numel(N_histbin));
    dispnames=cell(size(rmax_Psi_vals));
    for i_r=1:numel(rmax_Psi_vals)
        barmat(i_r,:)=reshape(N_hist(i_T,i_r,:),size(N_histbin));
        if rmax_Psi_vals(i_r)<0
            dispnames{i_r}="Delaunay";
        elseif rmax_Psi_vals(i_r)==0
            dispnames{i_r}="Topological";
        else
            dispnames{i_r}=sprintf("$r_{\\textrm{max}}=%.2f$",rmax_Psi_vals(i_r));
        end
            
    end
    
    figure
    plot(N_histbin,barmat);
    hleg{i_T}=legend(dispnames,'NumColumns',2,'Location','northwest','interpreter','latex');
    hleg{i_r}=legend('NumColumns',2,'Location','northeast','interpreter','latex');
    if i_T <= N_T
        title(sprintf("$N = %d,T = %.2f$",sqrtN,T_vals(i_T)),"interpreter",'latex');
    else
        title(sprintf("$N = %d$, %s",sqrtN,T_str(i_T)),"interpreter",'latex');
    end
    ylim([0 1.4]);
end


% for i_T = 1:N_T
%     figure
%     for i_r=3:numel(rmax_Psi_vals)
%         dispname=sprintf("$r_{\\textrm{max}}=%.2f$",rmax_Psi_vals(i_r));
%         bar(N_histbin,reshape(N_hist(i_T,i_r,:),size(N_histbin)),'grouped',...
%             'DisplayName',dispname);
%         hold on;
%     end
%     dispname="Delaunay";
%     bar(N_histbin,reshape(N_hist(i_T,1,:),size(N_histbin)),'grouped',...
%         'DisplayName',dispname);
%     dispname="Topological";
%     bar(N_histbin,reshape(N_hist(i_T,2,:),size(N_histbin)),'grouped',...
%         'DisplayName',dispname);
%     
%     hleg{i_r}=legend('NumColumns',2,'Location','northeast','interpreter','latex');
%     title(sprintf("$N = %d,T = %.2f$",sqrtN,T_vals(i_T)),"interpreter",'latex');
% end
% 
