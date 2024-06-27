function [nb,nb_n,dtr,r_new] = delaunay_nb(r,L,r_extend)
%DELAUNAY_NB Calculates neighbors via delaunay triangulation
%   Extends the prticle list by distance r_extend beyond the simulation box
%   of size L. Particle positions given by r.
%   Returns
%   nb     cell array with neighbors of the i-th particle
%   nb_n   array with number of neighbros 
%   Uses code from https://www.mathworks.com/matlabcentral/answers/163115-how-do-i-get-the-neighbors-of-a-vertex-in-a-delaunay-triangulation

nb=cell(1,numel(r)/2);
nb_n=zeros(1,numel(r)/2);

[r_new,~] = extend_pbc(r,[],L,r_extend);
x=r_new(1,:);
y=r_new(2,:);

% Calculates Delaunay triangulation
dtr=delaunayTriangulation(x(:),y(:));
% Calculates vertices attached to each particle.
tria_list=vertexAttachments(dtr);
for i = 1:numel(r)/2
    % 2. Use the connectivity list to get the vertex indices of all these
    % triangles
    vertices_cur = dtr.ConnectivityList(tria_list{i},:);
    % 3. Find all the unique vertices and remove the current vertex
    nb{i} = setdiff(unique(vertices_cur), i);
    nb_n(i)=numel(nb{i});
end

end

