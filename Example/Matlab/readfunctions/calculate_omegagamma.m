function calculate_omegagamma(infile_base) %, system, L)
% FUNCTION CALCULATE_OMEGAGAMMA calculates mean values and errors
% for the fitted values omega and gamma for the time correlation function
% C_mperpmperp(q,t) (gmperpmperp)

    %% Read data
    infile_std=sprintf('%s.mat',infile_base);
    infile_jkn=sprintf('%s_jackknife_gmperp.mat',infile_base);
    outfile=sprintf('%s_omegagamma_gmperp.mat',infile_base);
    load(infile_jkn,'qbin','t','jkn');
    load(infile_std,"gmperpmperp");
    if size(t,1) == 1
        t=t';
    end
    cf=jkn;
    n_t=numel(t);
    n_q=numel(qbin);
    n_runs=size(jkn,1);

    %% Initialize functions and vectors
    fitfunc="exp(-gamma*x/2)*(cos(omega_1*x)+gamma/omega_1/2*sin(omega_1*x))";
    fitobs=cell(n_runs,n_q);
    
    gamma_jkn=zeros(n_runs,n_q);
    omega_1_jkn=gamma_jkn;
    gamma_cimax_jkn=gamma_jkn;
    gamma_cimin_jkn=gamma_jkn;
    omega_1_cimax_jkn=gamma_jkn;
    omega_1_cimin_jkn=gamma_jkn;

    cf_mean=real(reshape(gmperpmperp,n_q,n_t));
    startpoint=[1e-4,1e-3];
    gamma_mean=zeros(1,n_q);
    omega_1_mean=0*gamma_mean;
    gamma_cimin_mean=0*gamma_mean;
    gamma_cimax_mean=0*gamma_mean;
    omega_1_cimin_mean=0*gamma_mean;
    omega_1_cimax_mean=0*gamma_mean;
    fitobs_mean=cell(size(gamma_mean));

    %% Calculating coefficients for the mean value, that is cf_mean
    for i_q = 1:n_q
        cf_cur=real(reshape(cf_mean(i_q,:),size(t)));
        cf_cur=cf_cur/cf_cur(1);

        if cf_cur(1) > 0
            fitobs_mean{i_q} = fit(t,cf_cur,fitfunc, ...
                'StartPoint',startpoint,...
                'Lower',[0,0]);
            
            gamma_mean(i_q)=fitobs_mean{i_q}.gamma;
            omega_1_mean(i_q)=fitobs_mean{i_q}.omega_1;

            ci=confint(fitobs_mean{i_q});
            gamma_cimin_mean(i_q)=ci(1,1);
            gamma_cimax_mean(i_q)=ci(2,1);
            omega_1_cimin_mean(i_q)=ci(1,2);
            omega_1_cimax_mean(i_q)=ci(2,2);

            startpoint=[gamma_mean(i_q),omega_1_mean(i_q)];
        end
    end

    %% Calculating coefficients for each of the jackknife samples
    for i_runs = 1:n_runs
        fprintf('%d \n',i_runs);
        for i_q = 1:n_q
            cf_cur=real(reshape(cf(i_runs,i_q,:),size(t)));
            cf_cur=cf_cur/cf_cur(1);
            
            startpoint=[gamma_mean(i_q),omega_1_mean(i_q)];
        
            if cf_cur(1) > 0
                fitobs{i_runs,i_q} = fit(t,cf_cur,fitfunc, ...
                    'StartPoint',startpoint,...
                    'Lower',[0,0]);
                
                gamma_jkn(i_runs,i_q)=fitobs{i_runs,i_q}.gamma;
                omega_1_jkn(i_runs,i_q)=fitobs{i_runs,i_q}.omega_1;
    
                ci=confint(fitobs{i_runs,i_q});
                gamma_cimin_jkn(i_runs,i_q)=ci(1,1);
                gamma_cimax_jkn(i_runs,i_q)=ci(2,1);
                omega_1_cimin_jkn(i_runs,i_q)=ci(1,2);
                omega_1_cimax_jkn(i_runs,i_q)=ci(2,2);

            end
        end
    end

    %% Performing jackknife estimates of all quantities of interest
    gamma_bias=(n_runs-1)*(mean(gamma_jkn) - gamma_mean);
    gamma_bcmean=gamma_mean-gamma_bias;
    gamma_pseudovals=n_runs*gamma_mean - (n_runs-1)*gamma_jkn;
    gamma_jackvar=1/(n_runs*(n_runs-1))*sum((gamma_pseudovals-gamma_bcmean).^2);
    gamma_ste=sqrt(gamma_jackvar);

    omega_1_bias=(n_runs-1)*(mean(omega_1_jkn) - omega_1_mean);
    omega_1_bcmean=omega_1_mean-omega_1_bias;
    omega_1_pseudovals=n_runs*omega_1_mean - (n_runs-1)*omega_1_jkn;
    omega_1_jackvar=1/(n_runs*(n_runs-1))*sum((omega_1_pseudovals-omega_1_bcmean).^2);
    omega_1_ste=sqrt(omega_1_jackvar);

    gamma_cimin_bias=(n_runs-1)*(mean(gamma_cimin_jkn) - gamma_cimin_mean);
    gamma_cimin_bcmean=gamma_cimin_mean-gamma_cimin_bias;
    gamma_cimin_pseudovals=n_runs*gamma_cimin_mean - (n_runs-1)*gamma_cimin_jkn;
    gamma_cimin_jackvar=1/(n_runs*(n_runs-1))*sum((gamma_cimin_pseudovals-gamma_cimin_bcmean).^2);
    gamma_cimin_ste=sqrt(gamma_cimin_jackvar);

    gamma_cimax_bias=(n_runs-1)*(mean(gamma_cimax_jkn) - gamma_cimax_mean);
    gamma_cimax_bcmean=gamma_cimax_mean-gamma_cimax_bias;
    gamma_cimax_pseudovals=n_runs*gamma_cimax_mean - (n_runs-1)*gamma_cimax_jkn;
    gamma_cimax_jackvar=1/(n_runs*(n_runs-1))*sum((gamma_cimax_pseudovals-gamma_cimax_bcmean).^2);
    gamma_cimax_ste=sqrt(gamma_cimax_jackvar);

    omega_1_cimin_bias=(n_runs-1)*(mean(omega_1_cimin_jkn) - omega_1_cimin_mean);
    omega_1_cimin_bcmean=omega_1_cimin_mean-omega_1_cimin_bias;
    omega_1_cimin_pseudovals=n_runs*omega_1_cimin_mean - (n_runs-1)*omega_1_cimin_jkn;
    omega_1_cimin_jackvar=1/(n_runs*(n_runs-1))*sum((omega_1_cimin_pseudovals-omega_1_cimin_bcmean).^2);
    omega_1_cimin_ste=sqrt(omega_1_cimin_jackvar);

    omega_1_cimax_bias=(n_runs-1)*(mean(omega_1_cimax_jkn) - omega_1_cimax_mean);
    omega_1_cimax_bcmean=omega_1_cimax_mean-omega_1_cimax_bias;
    omega_1_cimax_pseudovals=n_runs*omega_1_cimax_mean - (n_runs-1)*omega_1_cimax_jkn;
    omega_1_cimax_jackvar=1/(n_runs*(n_runs-1))*sum((omega_1_cimax_pseudovals-omega_1_cimax_bcmean).^2);
    omega_1_cimax_ste=sqrt(omega_1_cimax_jackvar);

%     bias=(n-1)*(jkn_mean - meanval);
%     bcmean=meanval -  bias;
%     pseudovals=n*meanval - (n-1)*jkn;
%     jackvar=1/(n*(n-1))*sum((pseudovals-bcmean).^2);
%     stderr=sqrt(jackvar);


    fprintf('Saving to %s\n',outfile);
    save(outfile,'qbin',...
        'gamma_mean','gamma_bias','gamma_bcmean','gamma_ste','gamma_jkn',...
        'omega_1_mean','omega_1_bias','omega_1_bcmean','omega_1_ste','omega_1_jkn',...
        'gamma_cimin_mean','gamma_cimin_bias','gamma_cimin_bcmean','gamma_cimin_ste','gamma_cimin_jkn',...
        'gamma_cimax_mean','gamma_cimax_bias','gamma_cimax_bcmean','gamma_cimax_ste','gamma_cimax_jkn',...
        'omega_1_cimin_mean','omega_1_cimin_bias','omega_1_cimin_bcmean','omega_1_cimin_ste','omega_1_cimin_jkn',...
        'omega_1_cimax_mean','omega_1_cimax_bias','omega_1_cimax_bcmean','omega_1_cimax_ste','omega_1_cimax_jkn')



%     cf_full=reshape(gmperpmperp_collect,n_runs,n_q,n_t);
%     [bcmean,bias,stderr,jkn]=jackknife_estimate(cf_full,@mean);
% 
%     save(outfile,...
%         't','n_t','qbin','n_q','n_runs',...
%         'jkn','bcmean','bias','stderr');

%     [bcmean,~,stderr,jkn] = jackknife_estimate(SCF_Spin_av_collect,@mean);
%     SCF_mean=bcmean;
%     SCF_err=stderr;
%     SCF_jkn=jkn;
% 
%     Upsilon=1/2/L^2*(-(H_x_collect + H_y_collect) - 1/kT*(I_x_collect.^2 + I_y_collect.^2));
%     [bcmean,~,stderr,jkn] = jackknife_estimate(Upsilon,@mean);
%     Upsilon_mean=bcmean;
%     Upsilon_err=stderr;
%     Upsilon_jkn=jkn;
% 
%     f_eta=@(x) kT/2/pi/mean(x);
%     [bcmean,~,stderr,jkn] = jackknife_estimate(Upsilon,f_eta);
%     eta_Upsilon_mean=bcmean;
%     eta_Upsilon_err=stderr;
%     eta_Upsilon_jkn=jkn;
% %     eta_Upsilon_mean=kT/2/pi/Upsilon_mean;
% %     eta_Upsilon_err=max(abs(kT/2/pi/(Upsilon_mean+Upsilon_err) - eta_Upsilon_mean),...
% %         abs(kT/2/pi/(Upsilon_mean-Upsilon_err) - eta_Upsilon_mean));
% 
%     save(outfile,...
%         'SCF_mean','SCF_err','SCF_jkn',...
%         'Upsilon_mean','Upsilon_err','Upsilon_jkn',...
%         'eta_Upsilon_mean','eta_Upsilon_err','eta_Upsilon_jkn');
end
