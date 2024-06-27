function calculate_gr_spin(infile, outfile, system, L, dr,n_part)
    if (~ isfile(outfile))
        [r,v,te,w] = mxy_snapshot_extract(infile,'rpto',system);
        x=r(1,:);
        y=r(2,:);
        vx=v(1,:);
        vy=v(2,:);
        
        verbose=false;
        blksize=1e3;
        
        N=numel(te);
        m=[sum(cos(te)),sum(sin(te))]/numel(te);
        mnorm=norm(m);
        te_m=atan2(m(2),m(1));
        te=atan2(sin(te-te_m),cos(te-te_m));

        s_par_2=sum(cos(te).^2)/N;
        s_par_4=sum(cos(te).^4)/N;
        s_perp_2=sum(sin(te).^2)/N;
        s_perp_4=sum(sin(te).^4)/N;

        sum_cos_2te=sum(cos(2*te))/N;
        sum_cos_4te=sum(cos(4*te))/N;
        sum_cos_6te=sum(cos(6*te))/N;
        sum_cos_8te=sum(cos(8*te))/N;
        
        [gr_hist, Cm, Cm_par, Cm_perp, Cm_cross, Cw, Cvx, Cvy, rbin] = ...
            twopointcorr_spin( x,y,te, vx,vy, w, dr, L,L, blksize,n_part,verbose);
        save(outfile, 'gr_hist', 'Cm', 'Cm_par', 'Cm_perp', ...
            'Cm_cross', 'Cw', 'Cvx', 'Cvy', 'rbin',...
            'm','mnorm','N','L','dr',...
            's_par_2','s_par_4','s_perp_2','s_perp_4',...
            'sum_cos_2te','sum_cos_4te','sum_cos_6te','sum_cos_8te');
        [gr,rbin,rw]=twopointcorr(r(1,:),r(2,:),.01,100,0);
        save(outfile,'gr','rbin','rw');
    else
%         [r,v,te,w] = mxy_snapshot_extract(infile,'rpto',system);
%         x=r(1,:);
%         y=r(2,:);
%         vx=v(1,:);
%         vy=v(2,:);
%         
%         verbose=false;
%         blksize=1e3;
%         
%         N=numel(te);
%         m=[sum(cos(te)),sum(sin(te))]/numel(te);
%         mnorm=norm(m);
%         te_m=atan2(m(2),m(1));
%         te=atan2(sin(te-te_m),cos(te-te_m));
% 
%         s_par_2=sum(cos(te).^2)/N;
%         s_par_4=sum(cos(te).^4)/N;
%         s_perp_2=sum(sin(te).^2)/N;
%         s_perp_4=sum(sin(te).^4)/N;
% 
%         sum_cos_2te=sum(cos(2*te))/N;
%         sum_cos_4te=sum(cos(4*te))/N;
%         sum_cos_6te=sum(cos(6*te))/N;
%         sum_cos_8te=sum(cos(8*te))/N;
%         
%         load(outfile,'gr_hist', 'Cm', 'Cm_par', 'Cm_perp', ...
%             'Cm_cross', 'Cw', 'Cvx', 'Cvy', 'rbin',...
%             'm','mnorm');
%         save(outfile, 'gr_hist', 'Cm', 'Cm_par', 'Cm_perp', ...
%             'Cm_cross', 'Cw', 'Cvx', 'Cvy', 'rbin',...
%             'm','mnorm','N','L','dr',...
%             's_par_2','s_par_4','s_perp_2','s_perp_4',...
%             'sum_cos_2te','sum_cos_4te','sum_cos_6te','sum_cos_8te');
        fprintf('Outfile already exists: %s\n',outfile);
    end
    
end
