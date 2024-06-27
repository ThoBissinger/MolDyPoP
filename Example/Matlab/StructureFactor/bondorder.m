function [psi_k,N_N] = bondorder(r,L,k,rmax)
%BONDORDER Calculates the bond order parameter of a configuration
%   Bond order psi6 defined  by equation (8) in Engel et al. 2013, PRE.
%   Variables
%   r         Positions of particles (2 x N vector)
%   L         System size (double or 2 x 1 vector)
%   k         k-fold bond order
%   rmax      radius for neighbor calculation (optional)

    N = numel(r)/2;
    if numel(L) == 2
        rho = N/L(1)/L(2);
    else
        rho = N/L^2;
    end

    if nargin <=3
        rmax = "Delaunay";
    end
    if isnumeric(rmax)
        
    end
    % Copying for pbc  
    if ~isnumeric(rmax) && ((rmax == "Delaunay") || (rmax == "delaunay"))
        [nb,N_N,dtr,r_new] = delaunay_nb(r,L,2);
        psi_k=zeros(size(r,2),1);
        for i=1:numel(nb)
            rdist=r_new(:,nb{i})-r_new(:,i);
            psi_k(i)=mean(exp(1i*k*atan2(rdist(2,:),rdist(1,:))));
        end
    else
        if rmax == 0
            K = k+1;
            if k <= 6
                [r_new,~] = extend_pbc(r,[],L,3);
            elseif k <= 20
                [r_new,~] = extend_pbc(r,[],L,4);
            else
                [r_new,~] = extend_pbc(r,[],L,min(L)/2);
            end
        else
            K=2*ceil( pi*(rmax)^2*rho - 1); % doulbe the expected number of particles
            [r_new,~] = extend_pbc(r,[],L,rmax);
        end
        N_new = numel(r)/2;
    
        % knnsearch extracts the K nearest neighbors of each particle, 
        % including the particle itself.
        [Ind_NN,D_NN]=knnsearch(r_new',r','K',K,'Distance','Euclidean');
        Ind_NN = Ind_NN(:,2:end); % first entry is the particle itself, ignore self-interaction
        D_NN = D_NN(:,2:end); % first entry is the particle itself, ignore self-interaction
    
    %    Ind_NN = Ind_NN(D_NN <= 1);
    %    D_NN = D_NN(D_NN <= 1);
        
        rdiffs = reshape(r_new(:,Ind_NN),[2,size(Ind_NN)])-r;
        rnormed = rdiffs./reshape(D_NN,1,N,K-1);
    
    %     psi_pre=exp(1i*k*atan2(reshape(rnormed(2,:,:),N,K-1),reshape(rnormed(1,:,:),N,K-1)));
        psi_pre=exp(1i*k*atan2(reshape(rnormed(2,:,:),N,K-1),reshape(rnormed(1,:,:),N,K-1)));
        if rmax == 0
            psi_k=mean(psi_pre,2);
            N_N = ones(size(psi_k))*k;
        else
            N_N = sum(D_NN<=rmax,2);
            psi_k = sum(psi_pre .* (D_NN<=rmax),2)./N_N;
        end
        psi_k(isnan(psi_k)) = 0;
    end
end

