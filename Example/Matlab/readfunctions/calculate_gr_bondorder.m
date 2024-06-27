function calculate_gr_bondorder(infile, outfile, system, L, dr,dPsi, n_part)
%     if (~ isfile(outfile))
        [r,~,~,~] = mxy_snapshot_extract(infile,'rpto',system);
        x=r(1,:);
        y=r(2,:);
        
        verbose=false;
        blksize=1e3;
        
        k_vals=[2,3,4,6];
        rmax_vals=[0,.66,.73,.8,1.08];

        [gr_hist, Psivals, Gcorrs, rbin] = twopointcorr_bondorder( x,y, k_vals,rmax_vals, dr, L, blksize,n_part,verbose);
        
        Psibin=-dPsi:dPsi:1+dPsi;
        Psi_hist=zeros(numel(k_vals),numel(rmax_vals),numel(Psibin)-1);
        for i_k=1:numel(k_vals)
            for i_rmax=1:numel(rmax_vals)
                aux_hist=histcounts( abs(Psivals(i_k,i_rmax,:)),Psibin');
                Psi_hist(i_k,i_rmax,:)=aux_hist/numel(x)/(Psibin(2)-Psibin(1));
            end
        end
        K = 21;
        [r_new,~] = extend_pbc(r,[],L,5);
        [Ind_NN,D_NN]=knnsearch(r_new',r','K',K+1,'Distance','Euclidean');
        rmax=3;
        dr=.01;
        r_hist_bin=0:dr:rmax;
        r_hist=zeros(K,numel(r_hist_bin)-1);
        for i_K = 1:K
            r_hist(i_K,:)=histcounts( D_NN(:,i_K+1),r_hist_bin') / size(D_NN,1) / dr;
        end
        r_hist_bin=r_hist_bin(1:end-1);

        Psimean=mean(Psivals,3);
        Psiabsmean=mean(abs(Psivals),3);
        absPsimean=abs(Psimean);
        save(outfile, 'gr_hist', 'Gcorrs', 'rbin', ...
            'Psibin','Psi_hist',...
            'Psimean','Psiabsmean','absPsimean',...
            'k_vals', 'rmax_vals',...
            'r_hist_bin','r_hist');
%     else
%         fprintf('Outfile already exists: %s\n',outfile);
%     end
    
end
