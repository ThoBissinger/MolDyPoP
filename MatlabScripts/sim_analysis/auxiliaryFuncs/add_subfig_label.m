function h_text=add_subfig_label(ax_h,str,corner,xscale,yscale,fontsize)
%ADD_SUBFIG_LABEL Adds label (a,b, etc) to axis
%   Adds label to a corner of the axis handle passed
%   INPUT
%   ax_h      axis handle (can be gca)
%   str       name string for the label
%   corner    corner in which the label should appear
%   fontsize  size of the font
%   xscale    x axis scale ("lin" or "log")
%   yscale    y axis scale ("lin" or "log")

    if nargin < 5
        fontsize = 7;
        fonttype = 'cmr';
    end
    if nargin < 3
        xscale = 'lin';
        yscale = 'lin';
    end
    xlim=get(ax_h,'XLim');
    xlog=log(xlim);
    ylim=get(ax_h,'YLim');
    ylog=log(ylim);
    if xscale == 'lin'
        xl = min(xlim) + .02*diff(xlim);
        xr = max(xlim) - .02*diff(xlim);
    else
        xl = exp(min(xlog) + diff(xlog)*.02);
        xr = exp(max(xlog) - diff(xlog)*.02);
    end
    if yscale == 'lin'
        yd = min(ylim) + .02*diff(ylim);
        yu = max(ylim) - .02*diff(ylim);
    else
        yd = exp(min(ylog) + diff(ylog)*.02);
        yu = exp(max(ylog) - diff(ylog)*.02);
    end
    
    if (corner == "ne" || corner == "northeast")
        loc = [xr,yu];
        hor_al = "right";
        vert_al = "top";
    elseif (corner == "nw" || corner == "northwest")
        loc = [xl,yu];
        hor_al = "left";
        vert_al = "top";
    elseif (corner == "sw" || corner == "southwest")
        loc = [xl,yd];
        hor_al = "left";
        vert_al = "bottom";
    elseif (corner == "se" || corner == "southeast")
        loc = [xr,yd];
        hor_al = "right";
        vert_al = "bottom";
    end
    h_text=text('Position',loc, 'String', str, ...
        'interpreter','latex',...
        'HorizontalAlignment',hor_al,...
        'VerticalAlignment',vert_al,...
        'Color','black','FontSize', fontsize,...
        'FontName','cmr12' );

end

