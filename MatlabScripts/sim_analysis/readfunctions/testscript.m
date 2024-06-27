infile_base='/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00/sqrtN_64/T_.14/samp_Dynamics';
infile_std=sprintf('%s.mat',infile_base);
infile_jkn=sprintf('%s_jackknife_gmperp.mat',infile_base);
infile_par=sprintf('%s_omegagamma_gmperp.mat',infile_base);
outfile=sprintf('%s_jackknife_gmperp.mat',infile_base);
load(infile_std,'gmperpmperp');
load(infile_jkn,'qbin','t','jkn');
S_par=load(infile_par);

n_t=numel(t);
n_q=numel(qbin);
n_runs=size(jkn,1);
res_factor=3/4*log(n_runs);

res_function=@(t,tau) exp(-t.^2/2./tau.^2);

cf_mean=reshape(real(gmperpmperp),n_q,n_t);
cf_mean=cf_mean./cf_mean(:,1);

%% Calculate FTs
tau_vals=res_factor./S_par.gamma_mean;
res_vals=res_function(reshape(t,n_t,1),reshape(tau_vals,1,n_q));

[ft_re_mean,om_vals]=FT_correlation(t, real(cf_mean)' .* res_vals , 0);
ft_re_mean=ft_re_mean.';
[ft_im_mean,~]=FT_correlation(t, imag(cf_mean)' .* res_vals , 0);
ft_im_mean=ft_im_mean.';

[ft_re_jkn,~]=FT_correlation(t, real(permute(jkn,[3 2 1]) .* res_vals) , 0);
ft_re_jkn=permute(ft_re_jkn,[3 2 1]);
[ft_im_jkn,~]=FT_correlation(t, imag(permute(jkn,[3 2 1]) .* res_vals) , 0);
ft_im_jkn=permute(ft_im_jkn,[3 2 1]);

ft_re_bias=(n_runs-1)*(reshape(mean(ft_re_jkn),n_q,n_t) - ft_re_mean);
ft_re_bcmean=ft_re_mean-ft_re_bias;
ft_re_pseudovals=n_runs*reshape(ft_re_mean,1,n_q,n_t) - (n_runs-1)*ft_re_jkn;
ft_re_jackvar=1/(n_runs*(n_runs-1))*sum((ft_re_pseudovals-reshape(ft_re_bcmean,1,n_q,n_t)).^2);
ft_re_ste=sqrt(ft_re_jackvar);

ft_im_bias=(n_runs-1)*(reshape(mean(ft_im_jkn),n_q,n_t) - ft_im_mean);
ft_im_bcmean=ft_im_mean-ft_im_bias;
ft_im_pseudovals=n_runs*reshape(ft_im_mean,1,n_q,n_t) - (n_runs-1)*ft_im_jkn;
ft_im_jackvar=1/(n_runs*(n_runs-1))*sum((ft_im_pseudovals-reshape(ft_im_bcmean,1,n_q,n_t)).^2);
ft_im_ste=sqrt(ft_im_jackvar);
%     cf_full=reshape(gmperpmperp_collect,n_runs,n_q,n_t);
%     [bcmean,bias,stderr,jkn]=jackknife_estimate(cf_full,@mean);


%% Initialize functions and vectors for fits
fitfunc="2*gamma*omega_0/((x^2 - omega_0^2)^2 + gamma * x^2)";
fitobs=cell(n_runs,n_q);

gamma_jkn=zeros(n_runs,n_q);
omega_0_jkn=gamma_jkn;
gamma_cimax_jkn=gamma_jkn;
gamma_cimin_jkn=gamma_jkn;
omega_0_cimax_jkn=gamma_jkn;
omega_0_cimin_jkn=gamma_jkn;

cf_mean=real(reshape(gmperpmperp,n_q,n_t));
startpoint=[1e-4,1e-3];
gamma_mean=zeros(1,n_q);
omega_0_mean=0*gamma_mean;
gamma_cimin_mean=0*gamma_mean;
gamma_cimax_mean=0*gamma_mean;
omega_0_cimin_mean=0*gamma_mean;
omega_0_cimax_mean=0*gamma_mean;
fitobs_mean=cell(size(gamma_mean));

ft_re_ompeak_mean=zeros(size(qbin));
ft_re_peak_mean=zeros(size(qbin));
ft_re_spline_peak_mean=zeros(size(qbin));
ft_re_spline_ompeak_mean=zeros(size(qbin));
ft_re_hwhm_mean=zeros(size(qbin));
ft_re_ompeak_jkn=zeros(n_runs,n_q);
ft_re_peak_jkn=zeros(n_runs,n_q);
ft_re_spline_peak_jkn=zeros(n_runs,n_q);
ft_re_spline_ompeak_jkn=zeros(n_runs,n_q);
ft_re_hwhm_jkn=zeros(n_runs,n_q);



%% Calculating coefficients for the mean value, that is ft_re_mean
for i_q = 1:n_q
    chi_cur=real(cf_mean(i_q,1));
    ft_cur=real(reshape(ft_re_mean(i_q,:),size(om_vals))) / chi_cur;
    if chi_cur > 0
        startpoint=[S_par.gamma_mean(i_q),S_par.omega_1_mean(i_q)];
        fitobs_mean{i_q} = fit(om_vals(:),ft_cur(:),fitfunc, ...
            'StartPoint',startpoint,...
            'Lower',[0,0]);
        
        gamma_mean(i_q)=fitobs_mean{i_q}.gamma;
        omega_0_mean(i_q)=fitobs_mean{i_q}.omega_0;

        ci=confint(fitobs_mean{i_q});
        gamma_cimin_mean(i_q)=ci(1,1);
        gamma_cimax_mean(i_q)=ci(2,1);
        omega_0_cimin_mean(i_q)=ci(1,2);
        omega_0_cimax_mean(i_q)=ci(2,2);

        [ft_max,i_max]=max(ft_cur);
        ft_re_ompeak_mean(i_q)=abs(om_vals(i_max));
        ft_re_peak_mean(i_q)=ft_max*chi_cur;

        ft_spline_neg=spline(om_vals(:),-ft_cur(:));
        [splinemax,splinemaxarg]=fnmin(ft_spline_neg);
        ft_re_spline_peak_mean(i_q)=abs(splinemax)*chi_cur;
        ft_re_spline_ompeak_mean(i_q)=abs(splinemaxarg);
        ft_spline_halfheight=spline(om_vals(:),ft_cur(:)-abs(splinemax)/2);
        ft_spline_zeros=fnzeros(ft_spline_halfheight);
        if isempty(ft_spline_zeros)
            ft_re_hwhm_mean(i_q)=Inf;
        else
            ft_re_hwhm_mean(i_q)=max(abs(ft_spline_zeros(1,:)))-abs(splinemaxarg);
        end

    end
end


%% Calculating coefficients for the jackknife distribution, that is ft_re_jkn
for i_runs = 1:n_runs
    fprintf('%d \n',i_runs);
    for i_q = 1:n_q
        chi_cur=real(jkn(i_runs,i_q,1));
        ft_cur=real(reshape(ft_re_jkn(i_runs,i_q,:),size(om_vals))) / chi_cur;
        if chi_cur > 0
            startpoint=[S_par.gamma_jkn(i_runs,i_q),S_par.omega_1_jkn(i_runs,i_q)];
            fitobs{i_runs,i_q} = fit(om_vals(:),ft_cur(:),fitfunc, ...
                'StartPoint',startpoint,...
                'Lower',[0,0]);
            
            gamma_jkn(i_runs,i_q)=fitobs{i_runs,i_q}.gamma;
            omega_0_jkn(i_runs,i_q)=fitobs{i_runs,i_q}.omega_0;
    
            ci=confint(fitobs_mean{i_q});
            gamma_cimin_jkn(i_runs,i_q)=ci(1,1);
            gamma_cimax_jkn(i_runs,i_q)=ci(2,1);
            omega_0_cimin_jkn(i_runs,i_q)=ci(1,2);
            omega_0_cimax_jkn(i_runs,i_q)=ci(2,2);
    
            [ft_max,i_max]=max(ft_cur);
            ft_re_ompeak_jkn(i_runs,i_q)=abs(om_vals(i_max));
            ft_re_peak_jkn(i_runs,i_q)=ft_max*chi_cur;
    
            ft_spline_neg=spline(om_vals(:),-ft_cur(:));
            [splinemax,splinemaxarg]=fnmin(ft_spline_neg);
            ft_re_spline_peak_jkn(i_runs,i_q)=abs(splinemax)*chi_cur;
            ft_re_spline_ompeak_jkn(i_runs,i_q)=abs(splinemaxarg);
            ft_spline_halfheight=spline(om_vals(:),ft_cur(:)-abs(splinemax)/2);
            ft_spline_zeros=fnzeros(ft_spline_halfheight);
            if isempty(ft_spline_zeros)
                ft_re_hwhm_jkn(i_runs,i_q)=Inf;
            else
                ft_re_hwhm_jkn(i_runs,i_q)=max(abs(ft_spline_zeros(1,:)))-abs(splinemaxarg);
            end
    
        end
    end
end
%%
ft_re_bias=(n_runs-1)*(reshape(mean(ft_re_jkn),n_q,n_t) - ft_re_mean);
  ft_re_bcmean=ft_re_mean-ft_re_bias;
  ft_re_pseudovals=n_runs*reshape(ft_re_mean,1,n_q,n_t) - (n_runs-1)*ft_re_jkn;
  ft_re_jackvar=1/(n_runs*(n_runs-1))*sum((ft_re_pseudovals-reshape(ft_re_bcmean,1,n_q,n_t)).^2);
  ft_re_ste=sqrt(ft_re_jackvar);

ft_im_bias=(n_runs-1)*(reshape(mean(ft_im_jkn),n_q,n_t) - ft_im_mean);
  ft_im_bcmean=ft_im_mean-ft_im_bias;
  ft_im_pseudovals=n_runs*reshape(ft_im_mean,1,n_q,n_t) - (n_runs-1)*ft_im_jkn;
  ft_im_jackvar=1/(n_runs*(n_runs-1))*sum((ft_im_pseudovals-reshape(ft_im_bcmean,1,n_q,n_t)).^2);
  ft_im_ste=sqrt(ft_im_jackvar);

  %
gamma_bias=(n_runs-1)*(reshape(mean(gamma_jkn),size(qbin)) - gamma_mean);
  gamma_bcmean=gamma_mean-gamma_bias;
  gamma_pseudovals=n_runs*reshape(gamma_mean,size(qbin)) - (n_runs-1)*gamma_jkn;
  gamma_jackvar=1/(n_runs*(n_runs-1))*sum((gamma_pseudovals-reshape(gamma_bcmean,size(qbin))).^2);
  gamma_ste=sqrt(gamma_jackvar);

omega_0_bias=(n_runs-1)*(reshape(mean(omega_0_jkn),size(qbin)) - omega_0_mean);
  omega_0_bcmean=omega_0_mean-omega_0_bias;
  omega_0_pseudovals=n_runs*reshape(omega_0_mean,size(qbin)) - (n_runs-1)*omega_0_jkn;
  omega_0_jackvar=1/(n_runs*(n_runs-1))*sum((omega_0_pseudovals-reshape(omega_0_bcmean,size(qbin))).^2);
  omega_0_ste=sqrt(omega_0_jackvar);

gamma_cimin_bias=(n_runs-1)*(reshape(mean(gamma_cimin_jkn),size(qbin)) - gamma_cimin_mean);
  gamma_cimin_bcmean=gamma_cimin_mean-gamma_cimin_bias;
  gamma_cimin_pseudovals=n_runs*reshape(gamma_cimin_mean,size(qbin)) - (n_runs-1)*gamma_cimin_jkn;
  gamma_cimin_jackvar=1/(n_runs*(n_runs-1))*sum((gamma_cimin_pseudovals-reshape(gamma_cimin_bcmean,size(qbin))).^2);
  gamma_cimin_ste=sqrt(gamma_cimin_jackvar);

gamma_cimax_bias=(n_runs-1)*(reshape(mean(gamma_cimax_jkn),size(qbin)) - gamma_cimax_mean);
  gamma_cimax_bcmean=gamma_cimax_mean-gamma_cimax_bias;
  gamma_cimax_pseudovals=n_runs*reshape(gamma_cimax_mean,size(qbin)) - (n_runs-1)*gamma_cimax_jkn;
  gamma_cimax_jackvar=1/(n_runs*(n_runs-1))*sum((gamma_cimax_pseudovals-reshape(gamma_cimax_bcmean,size(qbin))).^2);
  gamma_cimax_ste=sqrt(gamma_cimax_jackvar);

omega_0_cimin_bias=(n_runs-1)*(reshape(mean(omega_0_cimin_jkn),size(qbin)) - omega_0_cimin_mean);
  omega_0_cimin_bcmean=omega_0_cimin_mean-omega_0_cimin_bias;
  omega_0_cimin_pseudovals=n_runs*reshape(omega_0_cimin_mean,size(qbin)) - (n_runs-1)*omega_0_cimin_jkn;
  omega_0_cimin_jackvar=1/(n_runs*(n_runs-1))*sum((omega_0_cimin_pseudovals-reshape(omega_0_cimin_bcmean,size(qbin))).^2);
  omega_0_cimin_ste=sqrt(omega_0_cimin_jackvar);

omega_0_cimax_bias=(n_runs-1)*(reshape(mean(omega_0_cimax_jkn),size(qbin)) - omega_0_cimax_mean);
  omega_0_cimax_bcmean=omega_0_cimax_mean-omega_0_cimax_bias;
  omega_0_cimax_pseudovals=n_runs*reshape(omega_0_cimax_mean,size(qbin)) - (n_runs-1)*omega_0_cimax_jkn;
  omega_0_cimax_jackvar=1/(n_runs*(n_runs-1))*sum((omega_0_cimax_pseudovals-reshape(omega_0_cimax_bcmean,size(qbin))).^2);
  omega_0_cimax_ste=sqrt(omega_0_cimax_jackvar);

ft_re_ompeak_bias=(n_runs-1)*(reshape(mean(ft_re_ompeak_jkn),size(qbin)) - ft_re_ompeak_mean);
  ft_re_ompeak_bcmean=ft_re_ompeak_mean-ft_re_ompeak_bias;
  ft_re_ompeak_pseudovals=n_runs*reshape(ft_re_ompeak_mean,size(qbin)) - (n_runs-1)*ft_re_ompeak_jkn;
  ft_re_ompeak_jackvar=1/(n_runs*(n_runs-1))*sum((ft_re_ompeak_pseudovals-reshape(ft_re_ompeak_bcmean,size(qbin))).^2);
  ft_re_ompeak_ste=sqrt(ft_re_ompeak_jackvar);

ft_re_peak_bias=(n_runs-1)*(reshape(mean(ft_re_peak_jkn),size(qbin)) - ft_re_peak_mean);
  ft_re_peak_bcmean=ft_re_peak_mean-ft_re_peak_bias;
  ft_re_peak_pseudovals=n_runs*reshape(ft_re_peak_mean,size(qbin)) - (n_runs-1)*ft_re_peak_jkn;
  ft_re_peak_jackvar=1/(n_runs*(n_runs-1))*sum((ft_re_peak_pseudovals-reshape(ft_re_peak_bcmean,size(qbin))).^2);
  ft_re_peak_ste=sqrt(ft_re_peak_jackvar);

ft_re_spline_peak_bias=(n_runs-1)*(reshape(mean(ft_re_spline_peak_jkn),size(qbin)) - ft_re_spline_peak_mean);
  ft_re_spline_peak_bcmean=ft_re_spline_peak_mean-ft_re_spline_peak_bias;
  ft_re_spline_peak_pseudovals=n_runs*reshape(ft_re_spline_peak_mean,size(qbin)) - (n_runs-1)*ft_re_spline_peak_jkn;
  ft_re_spline_peak_jackvar=1/(n_runs*(n_runs-1))*sum((ft_re_spline_peak_pseudovals-reshape(ft_re_spline_peak_bcmean,size(qbin))).^2);
  ft_re_spline_peak_ste=sqrt(ft_re_spline_peak_jackvar);

ft_re_spline_ompeak_bias=(n_runs-1)*(reshape(mean(ft_re_spline_ompeak_jkn),size(qbin)) - ft_re_spline_ompeak_mean);
  ft_re_spline_ompeak_bcmean=ft_re_spline_ompeak_mean-ft_re_spline_ompeak_bias;
  ft_re_spline_ompeak_pseudovals=n_runs*reshape(ft_re_spline_ompeak_mean,size(qbin)) - (n_runs-1)*ft_re_spline_ompeak_jkn;
  ft_re_spline_ompeak_jackvar=1/(n_runs*(n_runs-1))*sum((ft_re_spline_ompeak_pseudovals-reshape(ft_re_spline_ompeak_bcmean,size(qbin))).^2);
  ft_re_spline_ompeak_ste=sqrt(ft_re_spline_ompeak_jackvar);

ft_re_hwhm_bias=(n_runs-1)*(reshape(mean(ft_re_hwhm_jkn),size(qbin)) - ft_re_hwhm_mean);
  ft_re_hwhm_bcmean=ft_re_hwhm_mean-ft_re_hwhm_bias;
  ft_re_hwhm_pseudovals=n_runs*reshape(ft_re_hwhm_mean,size(qbin)) - (n_runs-1)*ft_re_hwhm_jkn;
  ft_re_hwhm_jackvar=1/(n_runs*(n_runs-1))*sum((ft_re_hwhm_pseudovals-reshape(ft_re_hwhm_bcmean,size(qbin))).^2);
  ft_re_hwhm_ste=sqrt(ft_re_hwhm_jackvar);




save(outfile,...
    'om_vals','qbin',...
    'ft_re_mean','ft_re_bias','ft_re_bcmean','ft_re_ste','ft_re_jkn',...
    'ft_im_mean','ft_im_bias','ft_im_bcmean','ft_im_ste','ft_im_jkn',...
    'gamma_mean','gamma_bias','gamma_bcmean','gamma_ste','gamma_jkn',...
    'omega_0_mean','omega_0_bias','omega_0_bcmean','omega_0_ste','omega_0_jkn',...
    'gamma_cimin_mean','gamma_cimin_bias','gamma_cimin_bcmean','gamma_cimin_ste','gamma_cimin_jkn',...
    'gamma_cimax_mean','gamma_cimax_bias','gamma_cimax_bcmean','gamma_cimax_ste','gamma_cimax_jkn',...
    'omega_0_cimin_mean','omega_0_cimin_bias','omega_0_cimin_bcmean','omega_0_cimin_ste','omega_0_cimin_jkn',...
    'omega_0_cimax_mean','omega_0_cimax_bias','omega_0_cimax_bcmean','omega_0_cimax_ste','omega_0_cimax_jkn',...
    'ft_re_ompeak_mean','ft_re_ompeak_bias','ft_re_ompeak_bcmean','ft_re_ompeak_ste','ft_re_ompeak_jkn',...
    'ft_re_peak_mean','ft_re_peak_bias','ft_re_peak_bcmean','ft_re_peak_ste','ft_re_peak_jkn',...
    'ft_re_spline_peak_mean','ft_re_spline_peak_bias','ft_re_spline_peak_bcmean','ft_re_spline_peak_ste','ft_re_spline_peak_jkn',...
    'ft_re_spline_ompeak_mean','ft_re_spline_ompeak_bias','ft_re_spline_ompeak_bcmean','ft_re_spline_ompeak_ste','ft_re_spline_ompeak_jkn',...
    'ft_re_hwhm_mean','ft_re_hwhm_bias','ft_re_hwhm_bcmean','ft_re_hwhm_ste','ft_re_hwhm_jkn');