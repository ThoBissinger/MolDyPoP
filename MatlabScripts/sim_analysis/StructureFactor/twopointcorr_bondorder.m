%[gr_hist, Psivals, Gcorrs, rbin] = twopointcorr_bondorder( x,y,k_vals,rmax_vals, dr, L, blksize,n_part,verbose)
%
%   Computes the two point correlation function of a 2D lattice
%   of a fixed width and height, also correlations of the bond order
%   parameters
%
%   x - list of x coordinates of points
%   y - list of y coordinates of points
%   k_vals - values of k-fold symmetry
%   rmax_vals - maximum distance r for nearest neighbor calculations
%   dr - binning distance for histogram
%   L - (square) box size 
%   blksize - optional, default 1000, number of points to be considered
%   in a single step.
%   n_part - number of focal particles that are averaged over. Default is
%   all the particles
%   verbose - if true, will print which step is currently processed
%
%   gr_hist - histogram of gr
%   S - struct of bond order parameters
%   r - r-coordinates
%
%   Developed by Ilya Valmianski, alterations by Thomas Bissinger
%   For details on the definitions of the correlations, see thesis Thomas
%   Bissinger
%   email: ivalmian@ucsd.edu
function [gr_hist, Psivals, Gcorrs, rbin] = ...
    twopointcorr_bondorder( x,y, k_vals,rmax_vals, dr, L, blksize,n_part,verbose)
    
    

    %validate input     
    if length(x)~=length(y)
        error('Length of x should be same as length of y'); 
    elseif  numel(dr)~=1
        error('dr needs to have numel==1');
    elseif numel(x)~=length(x) || numel(y)~=length(y)
        error('Require x and y to be 1D arrays');
    end
    
    x = reshape(squeeze(x),[length(x) 1]);
    y = reshape(squeeze(y),[length(y) 1]);
    r = [x';y'];
    
    %validate/set blksize
    if nargin < 7
        blksize = 1000;
    elseif numel(blksize)~=1
        error('blksize must have numel = 1');
    elseif blksize < 1
        blksize = 1;
    elseif isinf(blksize) || isnan(blksize)
        blksize = length(x);
    end
    
    %validate/set n_part
    if nargin < 8
        n_part = numel(x);
    elseif n_part < 1 || n_part > numel(x)
        n_part = numel(x);
    end   
    %validate/set verbose
    if nargin ~= 9
        verbose = false;
    elseif numel(verbose)~=1
        error('verbose must have numel = 1');
    end   
        
        
    
    %number of particles
    totalPart = length(x);
    
    %largest radius possible
    maxR = L/2;
    
    %r bins and area bins
    rbin = dr:dr:maxR;
    if rbin(end)<maxR
        rbin=[rbin,maxR];
    end
    av_dens = totalPart/L^2;
    rareas = ((2*pi*rbin* dr)*av_dens);
    
    %preallocate space for corrfun/rw
    gr_hist = rbin*0;
    Gcorrs = zeros(numel(k_vals),numel(rmax_vals),numel(rbin));

    Psivals=zeros(numel(k_vals),numel(rmax_vals),numel(x));
    for i_k = 1 :numel(k_vals)
        for i_rmax = 1 :numel(rmax_vals)
            if rmax_vals(i_rmax) > 0
                Psivals(i_k,i_rmax,:) = bondorder(r,L,k_vals(i_k),rmax_vals(i_rmax));
            elseif rmax_vals(i_rmax) == 0
                Psivals(i_k,i_rmax,:) = bondorder(r,L,k_vals(i_k),"Delaunay");
            end
        end
    end
    
    %number of steps to be considered
    numsteps = ceil(n_part / blksize);

    %shuffle particles if n_part < total number of particles
    %loop will only run up to n_part, but of course we will keep the full
    %vector as for the sub-selection of focal particles, we still calculate
    %the full correlation for each focal particle
    if n_part < totalPart
        new_order=randperm(totalPart)';
        x = x(new_order);
        y = y(new_order);
        r = [x';y'];
        for i_k = 1 :numel(k_vals)
            for i_rmax = 1 :numel(rmax_vals)
                Psivals(i_k,i_rmax,:) = Psivals(i_k,i_rmax,new_order);
            end
        end
    end
    
    for j = 1:numsteps
 
        %loop through all particles and compute the correlation function
        indi = (j-1)*blksize+1;
        indf = min(n_part,j*blksize);
        
        if verbose
            disp(['Step ' num2str(j) ' of ' num2str(numsteps) '. ' ...
                'Analyzing points ' num2str(indi) ' to '  num2str(indf)...
                ' of total ' num2str(n_part)]);
        end
        
%         [gr_histArr, CmArr, Cm_parArr, Cm_perpArr, Cm_crossArr, CwArr, CvxArr, CvyArr] = ...
%             arrayfun(@ (xj,yj,tej,vxj,vyj,wj) onePartCorr_spins(xj,yj,tej,vxj,vyj,wj,x,y,te,vx,vy,w,rbin,m,Lx,Ly),...
%             x(indi:indf),y(indi:indf),te(indi:indf), ...
%             vx(indi:indf),vy(indi:indf),w(indi:indf), ...
%             'UniformOutput',false);
        [gr_histArr, GcorrsArr] = ...
            arrayfun(@ (j) onePartCorr_bondorder(j,x,y,Psivals,rbin,L),... arrayfun(@ (j) onePartCorr_spins(j,x,y,te,vx,vy,w,rbin,m,Lx,Ly),...
            (indi:indf)', ...
            'UniformOutput',false);
        
        gr_hist =  gr_hist + sum(cell2mat(gr_histArr),1);
        MGcorr=permute(reshape(cell2mat(GcorrsArr),numel(k_vals),numel(indi:indf),numel(rmax_vals),numel(rbin)),[2 1 3 4]);
        Gcorrs = Gcorrs + reshape(sum(MGcorr,1),numel(k_vals),numel(rmax_vals),numel(rbin));
      
%         
%         plot(r(gr_hist~=0),Cm(gr_hist~=0)./gr_hist(gr_hist~=0))
%         hold on;
%         legend show;
    end
   
    gr_hist = gr_hist / n_part ./ rareas;
    Gcorrs = Gcorrs / n_part  ./ reshape(rareas,1,1,[]);
    
    
end





% function [gr_hist, Cm, Cm_par, Cm_perp, Cm_cross, Cw, Cvx, Cvy] = ...
%     onePartCorr_spins(xj,yj,tej,vxj,vyj,wj,x,y,te,vx,vy,w,r,m,Lx,Ly)
function [gr_hist, Gcorr] = ...
    onePartCorr_bondorder(j,x,y,Psivals,r,L)
% onePartCorr_spins(j,x,y,te,vx,vy,w,r,m,Lx,Ly)
    
    %Save coordinates of focal particle, as the particle order will be
    %changed in the process
    xj=x(j);
    yj=y(j);
    Psij=Psivals(:,:,j);
    
    Psisize=size(Psivals);

    %compute radiuses in the (xj,yj) centered coordinate system (periodic
    %box)
    xdiff=x-xj;
    xdiff(xdiff<-L/2)=xdiff(xdiff<-L/2)+L;
    xdiff(xdiff>L/2)=xdiff(xdiff>L/2)-L;

    ydiff=y-yj;
    ydiff(ydiff<-L/2)=ydiff(ydiff<-L/2)+L;
    ydiff(ydiff>L/2)=ydiff(ydiff>L/2)-L;

    rho=hypot(xdiff,ydiff)';
%     rho=rho(logical(rho))'; % removes the 0

    %compute maximum unbiased rho (beyond that, there are finite-size
    %effects)
    maxRho = L/2;
    
    %truncate to highest unbiased rho
    Psivals=Psivals(:,:,rho<maxRho);
    rho=rho(rho<maxRho);
    

    % histcounts sorts particles into their respective bins. The result is
    % (proportional to) gr_hist
    [gr_hist,edges,bin]=histcounts( rho,[-inf r]');
    
    % accumarray sums over particles according to the binning found in the
    % above call to histcounts

    for i_k = 1:Psisize(1)
        for i_rmax = 1:Psisize(2)
            Psicur=reshape(Psivals(i_k,i_rmax,:),[],1);
            Gcorr(i_k,i_rmax,:) = accumarray(bin', ...
                Psicur, ...
                size(gr_hist'),@sum).';
        end
    end
    

    % The above calculation included a self part which is now removed by
    % hand
    gr_hist(1) = gr_hist(1) - 1;
    Gcorr(:,:,1) = Gcorr(:,:,1) - reshape(Psij,Psisize(1),Psisize(2),1);

    % Correct weighting with the jth particle (the particle the sum is
    % centered around)
    Gcorr = Gcorr .* reshape(conj(Psij),Psisize(1),Psisize(2),1);
    
    

    %     for i_rho = 1:numel(rho)
%         i_r = bin(i_rho); 
%         gr_hist(i_r) = gr_hist(i_r) + 1;
%         Cm(i_r) = Cm(i_r) + cos(te(i_rho) - tej);
%         Cm_par(i_r) = Cm_par(i_r) + cos(te(i_rho));
%         Cm_perp(i_r) = Cm_perp(i_r) + sin(te(i_rho));
% %         Cm_cross(i_r) = Cm_cross(i_r) + cos(te(i_r));
%     end
%     Cm_cross = Cm_par * mnorm;
%     Cm_par = Cm_par * cos(tej);
%     Cm_perp = Cm_perp * sin(tej);

end

