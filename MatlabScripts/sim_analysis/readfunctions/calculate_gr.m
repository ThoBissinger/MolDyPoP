function calculate_gr(infile, outfile, system)
    if (~ isfile(outfile))
        [r,~,~,~] = mxy_snapshot_extract(infile,'r',system);
        [gr,rbin,rw]=twopointcorr(r(1,:),r(2,:),.01,100,0);
        save(outfile,'gr','rbin','rw');
    else
        fprintf('Outfile already exists: %s',infile);
    end
    
end
