addpath '/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/StructureFactor'
addpath '/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/ScanningSnapshots'

clear
T=.01;
sqrtN=128;
L=74;
% L=74 * [1;sqrt(3)/2];
infile='/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_128/T_.01/run_4/output/snapshot_Dynamics_final.out';
[r,~,te,~] = mxy_snapshot_extract(infile,'rt','mxy','s');
rho=sqrtN^2/L^2;

%%
rmax=1;
[r_new,~] = extend_pbc(r,[],L,2);

x=r_new(1,:);
y=r_new(2,:);
%%
dtr=delaunayTriangulation(x(:),y(:));
nb=cell(numel(r)/2);

tria_list=vertexAttachments(dtr);
for i = 1:numel(r)/2
    % 2. Use the connectivity list to get the vertex indices of all these
    % triangles
    vertices_cur = dtr.ConnectivityList(tria_list{i},:);
    % 3. Find all the unique vertices and remove the current vertex
    nb{i} = setdiff(unique(vertices_cur), i);
end


%%
figure
triplot(dtr.ConnectivityList,x,y);
xlim([16.5,19]);
ylim([43,48]);
pbaspect([1 1 1]);

