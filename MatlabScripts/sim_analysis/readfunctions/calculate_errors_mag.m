function calculate_errors_mag(infile_base, outfile) %, system, L)
% FUNCTION CALCULATE_ERRORS_MAG calculates mean values and errors
% associated with diffferent magnetization-related properties. Mostly the
% magnetization itself, its variance and the binder cumulant.
    infile_std=sprintf('%s_collect.mat',infile_base);
    load(infile_std,'absM_av_collect',"M_2_av_collect","M_4_av_collect");

    [bcmean,~,stderr,jkn] = jackknife_estimate(absM_av_collect,@mean);
    m_mean=bcmean;
    m_err=stderr;
    m_jkn=jkn;

    var_vals=M_2_av_collect-m_mean.^2;
    [bcmean,~,stderr,jkn] = jackknife_estimate(var_vals,@mean);
    mvar_mean=bcmean;
    mvar_err=stderr;
    mvar_jkn=jkn;

    binder_vals=1 - 1/3 * M_4_av_collect/mean(M_2_av_collect)^2;
    [bcmean,~,stderr,jkn] = jackknife_estimate(binder_vals,@mean);
    binder_mean=bcmean;
    binder_err=stderr;
    binder_jkn=jkn;

%     [bcmean,~,stderr,jkn] = jackknife_estimate(absM_av_collect,@var);
%     mvar_mean=bcmean;
%     mvar_err=stderr;
%     mvar_jkn=jkn;
% 
%     f_binder=@(x) 1 - 1/3 * mean(x.^4)/mean(x.^2).^2;
%     [bcmean,~,stderr,jkn] = jackknife_estimate(absM_av_collect,f_binder);
%     binder_mean=bcmean;
%     binder_err=stderr;
%     binder_jkn=jkn;

    save(outfile,'m_mean','m_err','m_jkn',...
        'mvar_mean','mvar_err','mvar_jkn',...
        'binder_mean','binder_err','binder_jkn');

    
%     infile_helicity=sprintf('%s_helicity.mat',infile_base);
%     load(infile_std,"absM","SCF_Spin_av","gmperpmperp","temperature","H",...
%         "rbin","averaging_times","qbin","N");
%     load(infile_helicity,"H_x","H_y","I_x","I_y");
% 
%     T_mean=mean(temperature);
%     T_var=var(temperature);
%     absM_mean=mean(absM);
%     absM_var=var(absM);
%     M_2_mean=mean(absM.^2);
%     M_2_var=var(absM.^2);
%     M_4_mean=mean(absM.^4);
%     M_4_var=var(absM.^4);
% 
%     Upsilon=1/2/L^2*(H_x+H_y - 1/T_mean*(I_x^2+I_y^2));
% 
%     eta_magfit=-4*log(absM_mean)/log(2*N);
%     i_fit=find(( rbin > 1).* (rbin < .4*L));
%     fitob_SCF_pow=fit(rbin(i_fit)',SCF_Spin_av(i_fit)',"a*x^(-eta)",...
%         "StartPoint",[1,.25],"Lower",[0,0]);
%     a_SCF_pow=fitob_SCF_pow.a;
%     eta_SCF_pow=fitob_SCF_pow.eta;
%     fitob_SCF_exp=fit(rbin(i_fit)',SCF_Spin_av(i_fit)',"a*exp(-x/xi)",...
%         "StartPoint",[1,1],"Lower",[0,0]);
%     a_SCF_exp=fitob_SCF_exp.a;
%     xi_SCF_exp=fitob_SCF_exp.xi;
% 
%     t=averaging_times;
%     q_integers=find((qbin(1:end-q_off)/qbin(1)>round(qbin(1:end-q_off)/qbin(1))-1e-3) ...
%         .* (qbin(1:end-q_off)/qbin(1)<round(qbin(1:end-q_off)/qbin(1))+1e-3));
%     elseif qmode == "qfull"
%         q_integers=1:numel(qbin)-q_off;
%     end
%     q_integers=find(gmperpmperp(q_integers)~=0);
%     
%     n_q = numel(qbin);
%     n_q_int = numel(q_integers);
%     ga_vals=zeros(1,n_q_int);
%     om_vals=zeros(1,n_q_int);
%     
%     q_vals=qbin(q_integers);
%     for i_q_int = 1:numel(q_vals)
%         i_q = q_integers(i_q_int);
%         fprintf('q = %.3f, ',q_vals(i_q_int));
%         cf=real(gmperpmperp(i_q:numel(qbin):end))/real(gmperpmperp(i_q));
%         fitob_sym=fit(t(:),cf(:),fitfunc_symm,'StartPoint',[1e-2, 1e-2], 'Lower', [0, 0]);
%         ga_vals(i_q_int)=fitob_sym.gamma;
%         om_vals(i_q_int)=abs(fitob_sym.omega_1);
%     

    
%     %     if (~ isfile(outfile))
%         if system == "xy_s"
%             [r,~,te,~] = mxy_snapshot_extract(infile,'rt','xy','s');
%         elseif system == "xy"
%             [r,~,te,~] = mxy_snapshot_extract(infile,'rt','xy','t');
%         else
%             [r,~,te,~] = mxy_snapshot_extract(infile,'rt',system);
%         end
%         if (system == "xy" || system == "xy_s") % For a peridoic system, a small offset is needed for proper boundary handling
%             r=r+[.5;.5];
%             jfunc=@(r) 1;
%             ufunc=@(r) 0;
%         else
%             jfunc=@(r) (1-r).^2;
%             ufunc=@(r) 4*(1-r).^2;
%         end
% 
%         [H_x,H_y,I_x,I_y,H_s,H_r] = helicitymodulus(r,te,L,jfunc,ufunc);
%         
%         save(outfile, 'H_x','H_y','I_x','I_y','H_s','H_r');
%     
end
