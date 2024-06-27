function h_annot=add_subfig_label_box(ax_h,str,corner,fontsize,fonttype,width,height)
%ADD_SUBFIG_LABEL_BOX Adds label (a,b, etc) to axis
%   Adds label to a corner of the axis handle passed
%   INPUT
%   ax_h      axis handle (can be gca)
%   str       name string for the label
%   corner    position (specifies corner)
%   fontsize  size of the font
%   fonttype  type of the font
%   width     width of the box containing label
%   height    height of the box containing label

    if nargin < 6
       width=.03;
       height=.05;
    end
    if nargin < 4
        fontsize = 7;
        fonttype = 'cmr';
    end
       
    pos = get(ax_h,'position');
    l_lim = pos(1);             % left limit of current axis
    r_lim = pos(1)+pos(3);      % right limit of current axis
    b_lim = pos(2);             % bottom limit of current axis
    t_lim = pos(2)+pos(4);      % top limit of current axis
    if (corner == "ne" || corner == "northeast")
        dim = [r_lim-width t_lim-height width height];
    elseif (corner == "nw" || corner == "northwest")
        dim = [l_lim t_lim-height width height];
    elseif (corner == "sw" || corner == "southwest")
        dim = [l_lim b_lim width height];
    elseif (corner == "se" || corner == "southeast")
        dim = [r_lim-width b_lim width height];
    end
    h_annot=annotation('textbox',dim,'String',str,...
        'Color','white',...
        'LineWidth',.5,...
        'Units','normalized',...
        'FitBoxToText','off',...
        'FontName',fonttype,...
        'FontSize',fontsize);
        

end

