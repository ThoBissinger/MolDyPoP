clear
cd /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/DensityCoarseGraining
run('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/initialization_script.m')
pathbase=pwd;
model="mxy";
sqrtN = 256;
L=74/128*sqrtN;
T_str = ".14";
T = str2double(T_str);
runnr = 1;
% S=load('/data/scc/thobi/220829_EquilibrationCheck/anneal/aligned/mxy_3.00/sqrtN_128/T_.01/samp_eq_collect',"absM_collect","averaging_times");
infile = "/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_" + sqrtN ...
    + "/T_" + T_str + "/run_" + runnr + "/output/snapshot_Dynamics_final.out";

[r,~,te,~] = mxy_snapshot_extract(infile,'rt','mxy');

img_res=600;

mesh_size_vals = [.5, .8, 1, 1.5, 2, 4, 6, 10];
N_mesh = numel(mesh_size_vals);

m=[sum(cos(te)),sum(sin(te))]/numel(te);
te_0=atan2(m(2),m(1));
absM=vecnorm(m);

%% Snapshot of anomalous configuration
img_size=[0 0 2*columnwidth_cm 2*columnwidth_cm];
figure
set(gcf,'units','centimeters','OuterPosition',img_size);
single_snap_fig(infile,'mxy','spins',"s");
xlim([0 L]);
ylim([0 L]);
axis off
% axis tight
% hXLabel = xlabel('x','interpreter','latex', ...
%     'FontName', 'cmr12','FontSize', fontsize_ax_labels);
% hYLabel = ylabel('$y$','interpreter','latex', ...
%     'FontName', 'cmr12','FontSize', fontsize_ax_labels);
figname=sprintf('%s/spinsnaps/%s_spins_%d_%.3f_%d',pathbase,model,sqrtN,T,runnr);
fprintf('Creating figure %s\n',figname)
exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',img_res);



%% Coarse grained densities
% for i_mesh = 1:N_mesh
for i_mesh = 1:N_mesh
    mesh_size=mesh_size_vals(i_mesh);
    cell_area=mesh_size.^2;
    [rho,centers]=hist3(r',[round(L/mesh_size),round(L/mesh_size)]);
    rho = rho/cell_area;
    
    x_lims=0:mesh_size:L;
    [X,Y] = meshgrid(x_lims(x_lims<L),x_lims(x_lims<L));

    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    set(gcf,'units','centimeters','InnerPosition',[0 0 columnwidth_cm columnwidth_cm]);
    hplot=pcolor(X,Y,rho);
    colorbar;
    colormap autumn;
    clim([2 4])
    if mesh_size >= 4
        hplot.LineWidth=.001;
    else
        hplot.LineStyle='none';
    end
    hplot.FaceColor='interp';
    rho_var=var(reshape(rho,[1 numel(rho)]));
    title_str=sprintf("$\\Delta r/\\sigma = %.1f$, $\\sqrt{\\langle \\Delta\\rho^2\\rangle} = %.3f$",mesh_size,sqrt(rho_var));
    htitle=title(title_str,"interpreter","latex","FontSize",14,...
        "FontName","cmr12");
    axis off
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    set(gcf,'units','centimeters','InnerPosition',[0 0 columnwidth_cm columnwidth_cm]);
    
    figname=sprintf('%s/densityCG/%s_rho_%d_%.3f_%d_%.1f',pathbase,model,sqrtN,T,runnr,mesh_size);
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',img_res/2);
%     exportgraphics(gcf,sprintf('%s.pdf',figname));

end

%% 
ind_func=@(x,y,dist) [ceil(x/dist);ceil(y/dist)];
for i_mesh = 1:N_mesh
    mesh_size=mesh_size_vals(i_mesh);
%     dist=L/mesh_size;
    cell_area=mesh_size.^2;

    ind_mat=ind_func(r(1,:),r(2,:),mesh_size);
    i_max=max(max(ind_mat));

    part_cell=cell(i_max,i_max);
    te_cg=zeros(i_max,i_max);
    rho_cg=te_cg;
    mpar_cg=te_cg;
    mperp_cg=te_cg;
%     part_in_cell_func=@(i,j) find((ind_mat(1,:) == i) .* (ind_mat(2,:) == j));
%     part_in_cell_func=@(v) find((ind_mat(1,:) == v(1)) .* (ind_mat(2,:) == v(2)));
%     part_cell=arrayfun(part_in_cell,[1:i_max],[1:i_max]);
    for i =1:i_max
        for j = 1:i_max
            part_cell{i,j} = find((ind_mat(1,:) == i) .* (ind_mat(2,:) == j));
            rho_cg(i,j)=numel(part_cell{i,j});
            te_cg(i,j)=sum(te(part_cell{i,j})-te_0);
            mpar_cg(i,j)=sum(cos(te(part_cell{i,j})-te_0));
            mperp_cg(i,j)=sum(sin(te(part_cell{i,j})-te_0));
        end
    end
    te_cg=te_cg./rho_cg;
    te_cg(isnan(te_cg))=0;
    mpar_cg=mpar_cg./rho_cg;
    mpar_cg(isnan(mpar_cg)) = 0;
    mperp_cg=mperp_cg./rho_cg;
    mperp_cg(isnan(mperp_cg)) = 0;
    rho_cg=rho_cg/cell_area;
    absm_cg = sqrt(mperp_cg.^2 + mpar_cg.^2);

    te_mean(i_mesh)=mean(te_cg(rho_cg ~= 0 ));
    mpar_mean(i_mesh)=mean(mpar_cg(rho_cg ~= 0 ));
    mperp_mean(i_mesh)=mean(mperp_cg(rho_cg ~= 0 ));
    absm_mean(i_mesh)=mean(absm_cg(rho_cg ~= 0 ));
    
    te_std(i_mesh)=std(te_cg(rho_cg ~= 0 ));
    mpar_std(i_mesh)=std(mpar_cg(rho_cg ~= 0 ));
    mperp_std(i_mesh)=std(mperp_cg(rho_cg ~= 0 ));
    absm_std(i_mesh)=std(absm_cg(rho_cg ~= 0 ));
    cmapsize=60;
    cmap=hsv(cmapsize);
    rgbcmp(1:cmapsize,1)=[linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3)];
    rgbcmp(1:cmapsize,2)=[zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3)];
    rgbcmp(1:cmapsize,3)=[linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3)];

    x_lims=0:mesh_size:L;
    [X,Y] = meshgrid(x_lims(x_lims<L),x_lims(x_lims<L));

    vars={te_cg;mpar_cg;mperp_cg;absm_cg};
    leg_names=["\theta","m_{\parallel}","m_{\perp}","m"];
    var_names=["te","mpar","mperp","absm"];
    c_lims={[-pi,pi],[-1,1],[-1,1],[0,1]};
    c_maps={rgbcmp,turbo,turbo,turbo};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % te_cg
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:numel(vars)
        figure
        curvar=vars{i};
        hplot=pcolor(X,Y,curvar);
        colorbar;
        colormap(c_maps{i});
    %     cb=colorbar;
        clim(c_lims{i})
        if mesh_size >= 4
            hplot.LineWidth=.001;
            hplot.FaceColor='interp';
        else
            hplot.LineStyle='none';
        end
        nonsing_vals=curvar(rho_cg~=0);
        std_vals(i)=std(nonsing_vals);
        mean_vals(i)=mean(nonsing_vals);
        title_str=sprintf("$\\Delta r/\\sigma = %.1f$, $\\sqrt{\\langle \\Delta %s\\rangle} = %.3f$",mesh_size,leg_names{i},std_vals(i));
        htitle=title(title_str,"interpreter","latex","FontSize",14,...
            "FontName","cmr12");
        axis off
        set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
        set(gcf,'units','centimeters','InnerPosition',[0 0 columnwidth_cm columnwidth_cm]);
    
        figname=sprintf('%s/%sCG/%s_%s_%d_%.3f_%d_%.1f',pathbase,var_names(i),model,var_names(i),sqrtN,T,runnr,mesh_size);
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',img_res/2);
    end


end


%%
aux=cat(2,kron((1:i_max)',ones(1,i_max)'),kron(ones(1,i_max)',(1:i_max)'));
S_res=arrayfun(@ (i,j) find((ind_mat(1,:) == i) .* (ind_mat(2,:) == j)),...
            aux(:,1), aux(:,2), ...
            'UniformOutput',false);