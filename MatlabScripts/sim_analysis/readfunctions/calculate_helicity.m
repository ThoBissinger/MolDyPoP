function calculate_helicity(infile, outfile, system, L)
%     if (~ isfile(outfile))
        if system == "xy_s"
            [r,~,te,~] = mxy_snapshot_extract(infile,'rt','xy','s');
        elseif system == "xy"
            [r,~,te,~] = mxy_snapshot_extract(infile,'rt','xy','t');
        else
            [r,~,te,~] = mxy_snapshot_extract(infile,'rt',system);
        end
        if (system == "xy" || system == "xy_s") % For a peridoic system, a small offset is needed for proper boundary handling
            r=r+[.5;.5];
            jfunc=@(r) 1;
            ufunc=@(r) 0;
        else
            jfunc=@(r) (1-r).^2;
            ufunc=@(r) 4*(1-r).^2;
        end

        [H_x,H_y,I_x,I_y,H_s,H_r] = helicitymodulus(r,te,L,jfunc,ufunc);
        
        save(outfile, 'H_x','H_y','I_x','I_y','H_s','H_r');
    
end
