function single_snap_fig(snapfile,system,figswitch,lattice)
%SINGLE_SNAP_FIG Creates an image of the snapshot in the snapfile
%   INPUT:
%   snapfile    location of snapshot file
%   system      the system to be studied
%               values:     'xy'      static xy model
%                           'mxy'     mobile xy model
%                           'vm'      Vicsek model
%                           'fvm'     frozen Vicsek model
%   figswitch   which kind of figure
%               values:     'spins'   spins quiver plot
%                           'pos'     positions of particle
%   lattice     lattice required for static system.
%               values:     's'     square lattice
%                           't'     trigonal lattice
%   OUTPUT:
%   r           positions in R^(2 x N)
%   p           momenta in R^(2 x N)
%   theta       angles in R^(1 x N)
%   omega       spin momenta in R^(1 x N)
    if nargin < 4
        lattice ='s';
    end
    [r,~,theta,~] = mxy_snapshot_extract(snapfile,'tr',system,lattice);
    L_x = ceil(max(r(1,:)));
    L_y = ceil(max(r(2,:)));
    

    n_angle=60;
    angle_res = 2 * pi / n_angle;
    cmapsize = n_angle;
    wrap_theta_vals=wrapToPi(theta);
    sx = .75*cos(theta);
    sy = .75*sin(theta);
%     colormap(jet);
%     c_map = hsv(n_angle);
%     cmap=hsv(cmapsize);
%     cmapsize=60;
    
    cmap=hsv(cmapsize);
    rgbcmp(1:cmapsize,1)=[linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3)];
    rgbcmp(1:cmapsize,2)=[zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3)];
    rgbcmp(1:cmapsize,3)=[linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3)];

%     rgbcmp(1:cmapsize,1)=[linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3)];
%     rgbcmp(1:cmapsize,2)=[zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3)];
%     rgbcmp(1:cmapsize,3)=[linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3)];

% bins=linspace(0,2*pi,cmapsize);
% angles=mod(angles,2*pi);
% col=discretize(angles,bins);

    axis([0 L_x 0 L_y])
    pbaspect([1 1 1])
    if (strcmp(figswitch,'spins') || strcmp(figswitch,'spins_pos'))
        for j = 1 : n_angle
            angle = - pi + (j - 1) * angle_res;
            indices=find((angle < wrap_theta_vals).*(wrap_theta_vals < angle + angle_res));
            quiver(r(1,indices),r(2,indices),...
                sx(indices),sy(indices),0,...
                'color',rgbcmp(j,:));
            if (strcmp(figswitch,'spins_pos'))
                scatter(r(1,:),r(2,:),5,'filled',...
                    'MarkerEdgeColor','b','MarkerFaceColor','b');
            end
%             axis([0 L_x 0 L_y])
            hold on;
        end
        hold off;
    elseif (strcmp(figswitch,'pos'))
        scatter(r(1,:),r(2,:),20,'filled');
        
    end
    pbaspect([1 1 1])
    fprintf('Lmax = %f\nabsM = %f\n', max(max(r)), sqrt(mean(cos(theta))^2 + mean(sin(theta))^2));
end

