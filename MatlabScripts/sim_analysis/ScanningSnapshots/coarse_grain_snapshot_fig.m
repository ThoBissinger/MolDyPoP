function coarse_grain_snapshot_fig(snapfile,system,resolution,fig_num,figswitch)
%COARSE_GRAIN_SNAPSHOT_FIG Creates an image of the snapshot in the snapfile
%after coarse-graining
%   INPUT:
%   snapfile    location of snapshot file
%   system      the system to be studied
%               values:     'xy'      static xy model
%                           'mxy'     mobile xy model
%   resolution  coarse graining resolution (size of averaging box)
%   fignum      number of the figure that should be created
%   figswitch   which kind of figure
%               values:     'spins'   spins quiver plot
%                           'pos'     positions of particle
%   OUTPUT:
%   r           positions in R^(2 x N)
%   p           momenta in R^(2 x N)
%   theta       angles in R^(1 x N)
%   omega       spin momenta in R^(1 x N)

    [r,~,theta,~] = mxy_snapshot_extract(snapfile,'tr',system);
    L_x = ceil(max(r(1,:)));
    L_y = ceil(max(r(2,:)));
    

    n_angle=30;
    angle_res = 2 * pi / n_angle;
    cmapsize = n_angle;
    wrap_theta_vals=wrapToPi(theta);
    sx = cos(theta);
    sy = sin(theta);
    
%     %% TODO
%     if (mod(L_x,resolution) == 0)
%         corase_rx = 0 : resolution : L_x;
%     else
%         corase_rx = [0 : resolution : L_x,L_x];
%     end
%     if (mod(L_y,resolution) == 0)
%         corase_ry = 0 : resolution : L_y;
%     else
%         corase_ry = [0 : resolution : L_y,L_y];
%     end
% %     [coarse_r(1,:),coarse_r(2,:)] = meshgrid(corase_rx,corase_ry);
%     [N,Xedges,Yedges,binX,binY] = histcounts2(r(1,:),r(2,:),corase_rx,corase_ry);
    [BinFrequency,Xedges,Yedges,binX,binY] = histcounts2(r(1,:),r(2,:),'BinWidth',[resolution resolution]);
    
    av_spins_x=zeros(size(BinFrequency));
    av_spins_y=zeros(size(BinFrequency));
    av_theta=zeros(size(BinFrequency));
    for i = 1:numel(BinFrequency(:,1))
        for j = 1:numel(BinFrequency(1,:))
            av_spins_x(i,j) = mean(sx(find((binX == i) .* (binY == j))));
            av_spins_y(i,j) = mean(sy(find((binX == i) .* (binY == j))));
            av_theta(i,j) = atan2(av_spins_y(i,j),av_spins_x(i,j));
        end
    end
    [rx,ry]=meshgrid(Xedges(1:end-1),Yedges(1:end-1));

    
%     surf(rx,ry,av_spins_x);
    
    
    
    
    
    
    colormap(jet);
%     c_map = hsv(n_angle);
    cmap=hsv(cmapsize);
    rgbcmp(1:cmapsize,1)=[linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3)];
    rgbcmp(1:cmapsize,2)=[zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3)];
    rgbcmp(1:cmapsize,3)=[linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3)];

% bins=linspace(0,2*pi,cmapsize);
% angles=mod(angles,2*pi);
% col=discretize(angles,bins);

    figure(fig_num)
    if (strcmp(figswitch,'spins'))
%         surf(rx,ry,(av_spins_x.^2 + av_spins_y.^2)');
        quiver(rx,ry,.6*resolution*av_spins_x',.6*resolution*av_spins_y',0);
%         for j = 1 : n_angle
%             angle = - pi + (j - 1) * angle_res;
%             indices=find((angle < wrap_theta_vals).*(wrap_theta_vals < angle + angle_res));
%             quiver(r(1,indices),r(2,indices),...
%                 sx(indices),sy(indices),0,...
%                 'color',rgbcmp(j,:));
% 
%             axis([0 L_x 0 L_y])
%             hold on;
%         end
        hold off;
    elseif (strcmp(figswitch,'pos'))
        plot(r(1,:),r(2,:),'o');
        axis([0 L_x 0 L_y])
    end
    
    fprintf('Lmax = %f\nabsM = %f\n', max(max(r)), sqrt(mean(cos(theta))^2 + mean(sin(theta))^2));
end

