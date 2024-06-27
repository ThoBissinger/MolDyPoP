function [r_new,te_new] = extend_pbc(r,te,L,sigma)
%EXTEND_PBC Extends coordinate data set by duplicating particles within
%sigma of the boundary on the other side of the boundary (i.e. particles
%with x < sigma will be duplicated at x' = x + L)
%   Variables
%   r         Positions of particles (2 x N vector)
%   te        Angles of particles (2 x N vector)
%   L         System size (double)
%   sigma     boundary with

    N = numel(r)/2;
    if numel(L) == 1
        L = [L;L];
    end

    % Finding particles at boundary
    i_x_l = find(r(1,:) < sigma); % left
    i_x_r = find(r(1,:) > L(1)-sigma); % right
    i_y_b = find(r(2,:) < sigma); % bottom
    i_y_t = find(r(2,:) > L(2)-sigma); % top

    r_new = [r,r(:,i_x_l)+[L(1);0],r(:,i_x_r)+[-L(1);0],... left and right
        r(:,i_y_b)+[0;L(2)],r(:,i_y_t)+[0;-L(2)],... bottom and top
        r(:,intersect(i_x_l,i_y_b))+[L(1);L(2)],... bottom left
        r(:,intersect(i_x_r,i_y_b))+[-L(1);L(2)],... bottom right
        r(:,intersect(i_x_l,i_y_t))+[L(1);-L(2)],... top left
        r(:,intersect(i_x_r,i_y_t))+[-L(1);-L(2)]... top right
        ];

    if ~isempty(te)
        te_new = [te,te(i_x_l),te(i_x_r),... left and right
            te(i_y_b),te(i_y_t),... bottom and top
            te(intersect(i_x_l,i_y_b)),... bottom left
            te(intersect(i_x_r,i_y_b)),... bottom right
            te(intersect(i_x_l,i_y_t)),... top left
            te(intersect(i_x_r,i_y_t))... top right
            ];
    else
        te_new=[];
    end
   
end

