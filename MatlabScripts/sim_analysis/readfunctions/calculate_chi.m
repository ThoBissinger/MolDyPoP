function calculate_chi(infile, outfile, system,L,qmax)
    if (~ isfile(outfile))
        if strcmp(system,"xy")
            lattice = 't';
        elseif strcmp(system,"xy_s")
            system = "xy";
            lattice = 's';
        else
            lattice='';
        end
        [r,~,te,~] = mxy_snapshot_extract(infile,'rt',system,lattice);
        [q_vals,chi,chimpar,chimperp,chi_mfree,chimpar_mfree,StructFac] = Susceptibility(r,te,L,qmax);
        save(outfile,'q_vals','chi','chimpar','chimperp', ...
            'chi_mfree','chimpar_mfree','StructFac');
    else
        fprintf('Outfile already exists: %s\n',infile);
    end
    
end
