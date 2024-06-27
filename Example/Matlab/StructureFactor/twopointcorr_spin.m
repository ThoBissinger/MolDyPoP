%[gr_hist, Cm, Cm_par, Cm_perp, Cm_cross, r] = twopointcorr_spin( x,y, te, dr, Lx,Ly, blksize,verbose)
%
%   Computes the two point correlation function of a 2D lattice
%   of a fixed width and height, also spin correlations
%
%   x - list of x coordinates of points
%   y - list of y coordinates of points
%   te - spin angle variable (theta)
%   dr - binning distance for histogram
%   Lx - box size in x-direction
%   Ly - box size in y-direction
%   blksize - optional, default 1000, number of points to be considered
%   in a single step.
%   n_part - number of focal particles that are averaged over. Default is
%   all the particles
%   verbose - if true, will print which step is currently processed
%
%   gr_hist - histogram of gr
%   Cm - total spin-spin correlation
%   Cm_par - parallel component spin-spin correlation
%   Cm_perp - perendicular component spin-spin correlation
%   Cm_cross - total magnetization to spin correlation
%   r - r-coordinates
%
%   Developed by Ilya Valmianski, alterations by Thomas Bissinger
%   For details on the definitions of the correlations, see thesis Thomas
%   Bissinger
%   email: ivalmian@ucsd.edu
function [gr_hist, Cm, Cm_par, Cm_perp, Cm_cross, Cw, Cvx, Cvy, rbin] = ...
    twopointcorr_spin( x,y,te, vx,vy, w, dr, Lx,Ly, blksize,n_part,verbose)
    
    m = 1/numel(te)*[sum(cos(te)),sum(sin(te))];
    te_m = atan2(m(2),m(1));
%     mperp = [m(2),-m(1)];
%     mnorm = m / norm(m);
%     mperp = [m(2),-m(1)];

    %validate input     
    if length(x)~=length(y),
        error('Length of x should be same as length of y'); 
    elseif  numel(dr)~=1,
        error('dr needs to have numel==1');
    elseif numel(x)~=length(x) || numel(y)~=length(y),
        error('Require x and y to be 1D arrays');
    end
    
    x = reshape(squeeze(x),[length(x) 1]);
    y = reshape(squeeze(y),[length(y) 1]);
    
    te = reshape(squeeze(te),[length(te) 1]);
%     te = atan2(sin(te - te_m),cos(te - te_m));

    vx = reshape(squeeze(vx),[length(vx) 1]);
    vy = reshape(squeeze(vy),[length(vy) 1]);
    w = reshape(squeeze(w),[length(w) 1]);
    
    %validate/set blksize
    if nargin < 7
        blksize = 1000;
    elseif numel(blksize)~=1,
        error('blksize must have numel = 1');
    elseif blksize < 1,
        blksize = 1;
    elseif isinf(blksize) || isnan(blksize),
        blksize = length(x);
    end
    
    %validate/set n_part
    if nargin < 8,
        n_part = numel(te);
    elseif n_part < 1 || n_part > numel(te)
        n_part = numel(te);
    end   
    %validate/set verbose
    if nargin ~= 9,
        verbose = false;
    elseif numel(verbose)~=1,
        error('verbose must have numel = 1');
    end   
        
        
    %real height/width
%     width = max(x)-min(x);
%     height = max(y)-min(y);
    
    %number of particles
    totalPart = length(x);
    
    %largest radius possible
    maxR = min(Lx/2,Ly/2);
    
    %r bins and area bins
    rbin = dr:dr:maxR;
    if rbin(end)<maxR;
        rbin=[rbin,maxR];
    end
    av_dens = totalPart/Lx/Ly;
    rareas = ((2*pi*rbin* dr)*av_dens);
    
    %preallocate space for corrfun/rw
    gr_hist = rbin*0;
    Cm = rbin*0;
    Cm_par = rbin*0;
    Cm_perp = rbin*0;
    Cm_cross = rbin*0;
    Cw = rbin*0;
    Cvx = rbin*0;
    Cvy = rbin*0;
    
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
        te = te(new_order);
        vx = vx(new_order);
        vy = vy(new_order);
        w = w(new_order);
    end
    
    for j = 1:numsteps
 
        %loop through all particles and compute the correlation function
        indi = (j-1)*blksize+1;
        indf = min(n_part,j*blksize);
        
        if verbose,
            disp(['Step ' num2str(j) ' of ' num2str(numsteps) '. ' ...
                'Analyzing points ' num2str(indi) ' to '  num2str(indf)...
                ' of total ' num2str(n_part)]);
        end
        
%         [gr_histArr, CmArr, Cm_parArr, Cm_perpArr, Cm_crossArr, CwArr, CvxArr, CvyArr] = ...
%             arrayfun(@ (xj,yj,tej,vxj,vyj,wj) onePartCorr_spins(xj,yj,tej,vxj,vyj,wj,x,y,te,vx,vy,w,rbin,m,Lx,Ly),...
%             x(indi:indf),y(indi:indf),te(indi:indf), ...
%             vx(indi:indf),vy(indi:indf),w(indi:indf), ...
%             'UniformOutput',false);
        [gr_histArr, CmArr, Cm_parArr, Cm_perpArr, Cm_crossArr, CwArr, CvxArr, CvyArr] = ...
            arrayfun(@ (j) onePartCorr_spins(j,x,y,te,vx,vy,w,rbin,m,Lx,Ly),...
            (indi:indf)', ...
            'UniformOutput',false);
        
        gr_hist =  gr_hist + sum(cell2mat(gr_histArr),1);
        Cm = Cm + sum(cell2mat(CmArr),1);
        Cm_par = Cm_par + sum(cell2mat(Cm_parArr),1);
        Cm_perp = Cm_perp + sum(cell2mat(Cm_perpArr),1);
        Cm_cross = Cm_cross + sum(cell2mat(Cm_crossArr),1);
        Cw = Cw + sum(cell2mat(CwArr),1);
        Cvx = Cvx + sum(cell2mat(CvxArr),1);
        Cvy = Cvy + sum(cell2mat(CvyArr),1);
      
%         
%         plot(r(gr_hist~=0),Cm(gr_hist~=0)./gr_hist(gr_hist~=0))
%         hold on;
%         legend show;
    end
   
    gr_hist = gr_hist / n_part ./ rareas;
    Cm = Cm / n_part  ./ rareas;
    Cm_par = Cm_par / n_part  ./ rareas;
    Cm_perp = Cm_perp / n_part  ./ rareas;
    Cm_cross = Cm_cross / n_part  ./ rareas;
    Cw = Cw / n_part  ./ rareas;
    Cvx = Cvx / n_part  ./ rareas;
    Cvy = Cvy / n_part  ./ rareas;
    
    
end





% function [gr_hist, Cm, Cm_par, Cm_perp, Cm_cross, Cw, Cvx, Cvy] = ...
%     onePartCorr_spins(xj,yj,tej,vxj,vyj,wj,x,y,te,vx,vy,w,r,m,Lx,Ly)
function [gr_hist, Cm, Cm_par, Cm_perp, Cm_cross, Cw, Cvx, Cvy] = ...
    onePartCorr_spins(j,x,y,te,vx,vy,w,r,m,Lx,Ly)
    
    %Save coordinates of focal particle, as the particle order will be
    %changed in the process
    xj=x(j);
    yj=y(j);
    tej=te(j);
    vxj=vx(j);
    vyj=vy(j);
    wj=w(j);
    %compute unit vectors
    mnorm = norm(m);

    %compute radiuses in the (xj,yj) centered coordinate system (periodic
    %box)
    xdiff=x-xj;
    xdiff(xdiff<-Lx/2)=xdiff(xdiff<-Lx/2)+Lx;
    xdiff(xdiff>Lx/2)=xdiff(xdiff>Lx/2)-Lx;

    ydiff=y-yj;
    ydiff(ydiff<-Ly/2)=ydiff(ydiff<-Ly/2)+Ly;
    ydiff(ydiff>Ly/2)=ydiff(ydiff>Ly/2)-Ly;

    rho=hypot(xdiff,ydiff)';
%     rho=rho(logical(rho))'; % removes the 0

    %compute maximum unbiased rho (beyond that, there are finite-size
    %effects)
    maxRho = min(Lx/2,Ly/2);
    
    %truncate to highest unbiased rho
    te=te(rho<maxRho);
    vx=vx(rho<maxRho);
    vy=vy(rho<maxRho);
    w=w(rho<maxRho);
    rho=rho(rho<maxRho);
    

    % histcounts sorts particles into their respective bins. The result is
    % (proportional to) gr_hist
    [gr_hist,edges,bin]=histcounts( rho,[-inf r]');
    
    % accumarray sums over particles according to the binning found in the
    % above call to histcounts
    Cm = accumarray(bin', ...
        cos(te-tej), ...
        size(gr_hist'),@sum)';
    Cm_par = accumarray(bin', ...
        cos(te), ...
        size(gr_hist'),@sum)';
    Cm_perp = accumarray(bin', ...
        sin(te), ...
        size(gr_hist'),@sum)';
    Cw = accumarray(bin', ...
        w, ...
        size(gr_hist'),@sum)';
    Cvx = accumarray(bin', ...
        vx, ...
        size(gr_hist'),@sum)';
    Cvy = accumarray(bin', ...
        vy, ...
        size(gr_hist'),@sum)';

    % The above calculation included a self part which is now removed by
    % hand
    gr_hist(1) = gr_hist(1) - 1;
    Cm(1) = Cm(1) - 1;
    Cm_par(1) = Cm_par(1) - cos(tej);
    Cm_perp(1) = Cm_perp(1) - sin(tej);
    Cw(1) = Cw(1) - wj;
    Cvx(1) = Cvx(1) - vxj;
    Cvy(1) = Cvy(1) - vyj;

    % Correct weighting with the jth particle (the particle the sum is
    % centered around)
    Cm_cross = Cm_par * mnorm;
    Cm_par = Cm_par * cos(tej);
    Cm_perp = Cm_perp * sin(tej);
    Cw = Cw * wj;
    Cvx = Cvx * vxj;
    Cvy = Cvy * vyj;

    

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

