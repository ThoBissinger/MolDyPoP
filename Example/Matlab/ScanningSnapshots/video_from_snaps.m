function video_from_snaps(snapfiles,times,system,figswitch)
%SINGLE_SNAP_FIG Creates an image of the snapshot in the snapfile
%   INPUT:
%   snapfiles   struct with all the files in the figure in R^(m)
%   times       double vector with all times in the figure in R^(m)
%   system      the system to be studied
%               values:     'xy'      static xy model
%                           'mxy'     mobile xy model
%   figswitch   which kind of figure
%               values:     'spins'   spins quiver plot
%                           'pos'     positions of particle

%     [r,~,theta,~] = mxy_snapshot_extract(snapfiles(1),'tr',system);
%     L_x = ceil(max(r(1,:)));
%     L_y = ceil(max(r(2,:)));

% t_vals = 1:999;
% for i = 1 : numel(t_vals)
%     fignames{i} = sprintf('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/src_eclipse/output/scaled_xy_%d.out', t_vals(i));
% end

    nsnap = numel(snapfiles);
    read_switch = 1;
    [r_initial,~,theta_initial,~] = mxy_snapshot_extract(snapfiles{1},'rt',system);
    N = numel(theta_initial);
    L_x = ceil(max(r_initial(1,:)));
    L_y = ceil(max(r_initial(2,:)));
    if (read_switch)
        if (strcmp(system,"mxy") || strcmp(system,"vm"))
            r_vals=zeros(nsnap,2,N);
        else
            r_vals = r_initial;
        end
        theta_vals=zeros(nsnap,N);
        for i = 1:nsnap
            disp(snapfiles{i});
            if (strcmp(system,"mxy") || strcmp(system,"vm"))
                [r_vals(i,:,:),~,theta_vals(i,:),~] = mxy_snapshot_extract(snapfiles{i},'tr',system);
            else
                [~,~,theta_vals(i,:),~] = mxy_snapshot_extract(snapfiles{i},'t',system);
            end
        end
        sx = .75*cos(theta_vals);
        sy = .75*sin(theta_vals);
        save('snap_data')
    else
        load('snap_data','r_vals', 'theta_vals', 'sx', 'sy')
    end

    
    n_angle=30;
    angle_res = 2 * pi / n_angle;
    cmapsize = n_angle;
    wrap_theta_vals=wrapToPi(theta_vals);
    
    colormap(jet);
%     c_map = hsv(n_angle);
%     cmap=hsv(cmapsize);
    rgbcmp(1:cmapsize,1)=[linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3)];
    rgbcmp(1:cmapsize,2)=[zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3)];
    rgbcmp(1:cmapsize,3)=[linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3)];

% bins=linspace(0,2*pi,cmapsize);
% angles=mod(angles,2*pi);
% col=discretize(angles,bins);

    figure(1)
    for i = 1 : nsnap
        if (strcmp(system,"mxy") || strcmp(system,"vm"))
            r = squeeze(r_vals(i,:,:));
        elseif ( i == 1)
            r = r_vals;
        end
        if (strcmp(figswitch,'spins'))
            for j = 1 : n_angle
                angle = - pi + (j - 1) * angle_res;
                indices=find((angle < wrap_theta_vals(i,:)).*(wrap_theta_vals(i,:) < angle + angle_res));
                
                
                quiver(r(1,indices),r(2,indices),...
                    sx(i,indices),sy(i,indices),0,...
                    'color',rgbcmp(j,:));

                axis([0 L_x 0 L_y])
                hold on;
            end
            hold off;
        elseif (strcmp(figswitch,'pos'))
            plot(r(1,:),r(2,:),'.');
            axis([0 L_x 0 L_y])
        end
        t = times(i);
        title(sprintf('NVE Simulation of %s model at t = %.2f', system, t));
        F(i) = getframe(gcf) ;
        drawnow
    end
    % create the video writer with 1 fps
    writerObj = VideoWriter('snapvid.avi');
    writerObj.FrameRate = 10;
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(F)
        % convert the image to a frame
        frame = F(i) ;    
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
    
    
%     fprintf('Lmax = %f\nabsM = %f\n', max(max(r)), sqrt(mean(cos(theta_vals(1,:)))^2 + mean(sin(theta_vals(1,:)))^2));
end

