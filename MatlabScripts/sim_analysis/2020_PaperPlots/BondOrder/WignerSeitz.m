sqrtN=16;
N=sqrtN^2;
T_str=".07";
infile="/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_" + sqrtN + "/T_" + T_str + "/run_4/output/snapshot_Dynamics_final.out";
system='mxy';
L=sqrtN/16*9.25;
K=20;
dr=.05;
[r,~,~,~]=mxy_snapshot_extract(infile,'r','mxy','t');
[r_new,~] = extend_pbc(r,[],L,2);

[vv,vn] = voronoin(r_new');
[Ind_NN,D_NN]=knnsearch(r_new',r','K',K,'Distance','Euclidean');

dtr=delaunayTriangulation(r_new(1,:)',r_new(2,:)');
%%
nb = cell(size(r,2),1);
for i=1:N
    nb_vec=[];
    for j=2:K
        s_int=numel(intersect(vn{i},vn{Ind_NN(i,j)}));
        if s_int > 1
            nb_vec = [nb_vec Ind_NN(i,j)];
        end
    end
    nb{i} = nb_vec;
end

%%
nb_dela=cell(numel(r)/2,1);
tria_list=vertexAttachments(dtr);
for i = 1:numel(r)/2
    % 2. Use the connectivity list to get the vertex indices of all these
    % triangles
    vertices_cur = dtr.ConnectivityList(tria_list{i},:);
    % 3. Find all the unique vertices and remove the current vertex
    nb_dela{i} = setdiff(unique(vertices_cur), i);
end


%%
close all;
i=1;
triplot(dtr);
hold on;

plot(r(1,i),r(2,i),'or',...
    'MarkerSize',20);
plot(r(1,nb{i}),r(2,nb{i}),'^b',...
    'MarkerSize',20);
xlim(r(1,i)+[-1,1]);
ylim(r(2,i)+[-1,1]);
pbaspect([ 1 1 1])


%%
for i = 1:N
    n1=reshape(sort(nb{i}),1,[]);
    n2=reshape(sort(nb_dela{i}),1,[]);
    ndiff(i)=max(abs(n1-n2));

end