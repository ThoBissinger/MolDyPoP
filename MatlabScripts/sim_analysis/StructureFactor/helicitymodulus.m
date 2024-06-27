function [H_x,H_y,I_x,I_y,H_s,H_r] = helicitymodulus(r,te,L,jfunc,ufunc)
%HELICITYMODULUS Calculates the helicity modulus of a configuration
%   Helicity modulus defined by equation (5) in Canova Levin Arenzon 2014.
%   Variables
%   r         Positions of particles (2 x N vector)
%   te        Angles of particles (2 x N vector)
%   be        Inverse temperature beta (double)
%   L         System size (double)
%   jfunc     function handle for the j interaction energy
%   ufunc     function handle for u interaction energy

    N = numel(te);
    K = 20;

    % Copying for pbc
    [r_new,te_new] = extend_pbc(r,te,L,1);
    N_new = numel(te);

    % knnsearch extracts the 15 nearest neighbors of each particle, 
    % including the particle itself. 15 was chosen such that (typically)
    % all nearest neighbors and not too much more are in the vector
    [Ind_NN,D_NN]=knnsearch(r_new',r','K',K,'Distance','Euclidean');
    Ind_NN = Ind_NN(:,2:end); % first entry is the particle itself, ignore self-interaction
    D_NN = D_NN(:,2:end); % first entry is the particle itself, ignore self-interaction

%    Ind_NN = Ind_NN(D_NN <= 1);
%    D_NN = D_NN(D_NN <= 1);
    
    rdiffs = reshape(r_new(:,Ind_NN),[2,size(Ind_NN)])-r;
        % easier computation. Removing D_NN > 1 entries first leads to a
        % lot of indexing issues
    tediffs = reshape(te_new(Ind_NN),size(Ind_NN))-te';
    
    indices=find(D_NN <= 1); % consider indices that agree with cutoff
    rdiffs=rdiffs(:,indices);
    tediffs=tediffs(indices)';
%     D_NN = D_NN(indices);
%     r_u = rdiffs./sqrt(rdiffs(1,:).^2 + rdiffs(2,:).^2); % unit vectors

%     origin_mat=kron((1:N_new),ones(1,K));
%     indices=find(D_NN <= 1); % ignore self-interaction and only consider particles within cutoff
%     rdiffs=r(:,Ind_NN(indices))-r(:,origin_mat(indices));
%     thetadiffs=te(Ind_NN(indices))-te(origin_mat(indices));
%     NN_x_component=rdiffs(1,:)./sqrt(rdiffs(1,:).^2 + rdiffs(2,:).^2);
%     NN_y_component=rdiffs(2,:)./sqrt(rdiffs(1,:).^2 + rdiffs(2,:).^2);
%    cosvals=cos(thetadiffs).*rdiffs;
    H_x = .5*sum(jfunc(vecnorm(rdiffs)).*cos(tediffs).*rdiffs(1,:).^2);
    H_y = .5*sum(jfunc(vecnorm(rdiffs)).*cos(tediffs).*rdiffs(2,:).^2);
    I_x = .5*sum(jfunc(vecnorm(rdiffs)).*sin(tediffs).*rdiffs(1,:));
    I_y = .5*sum(jfunc(vecnorm(rdiffs)).*sin(tediffs).*rdiffs(2,:));
    H_s = -.5*sum(jfunc(vecnorm(rdiffs)).*cos(tediffs));
    H_r = .5*sum(ufunc(vecnorm(rdiffs)));
   
end

