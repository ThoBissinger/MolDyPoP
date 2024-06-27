function calculate_errors_etafit(infile_base, outfile,L) %, system, L)
% FUNCTION CALCULATE_ERRORS_ETAFIT calculates mean values and errors
% with eta fits (for fits to static quantities)
    infile_std=sprintf('%s_collect.mat',infile_base);
    infile_helicity=sprintf('%s_helicity.mat',infile_base);
    load(infile_std,'SCF_Spin_av_collect','rbin');
    load(infile_helicity);

    [bcmean,~,stderr,jkn] = jackknife_estimate(SCF_Spin_av_collect,@mean);
    SCF_mean=bcmean;
    SCF_err=stderr;
    SCF_jkn=jkn;

    Upsilon=1/2/L^2*(-(H_x_collect + H_y_collect) - 1/kT*(I_x_collect.^2 + I_y_collect.^2));
    [bcmean,~,stderr,jkn] = jackknife_estimate(Upsilon,@mean);
    Upsilon_mean=bcmean;
    Upsilon_err=stderr;
    Upsilon_jkn=jkn;

    f_eta=@(x) kT/2/pi/mean(x);
    [bcmean,~,stderr,jkn] = jackknife_estimate(Upsilon,f_eta);
    eta_Upsilon_mean=bcmean;
    eta_Upsilon_err=stderr;
    eta_Upsilon_jkn=jkn;
%     eta_Upsilon_mean=kT/2/pi/Upsilon_mean;
%     eta_Upsilon_err=max(abs(kT/2/pi/(Upsilon_mean+Upsilon_err) - eta_Upsilon_mean),...
%         abs(kT/2/pi/(Upsilon_mean-Upsilon_err) - eta_Upsilon_mean));

    save(outfile,...
        'SCF_mean','SCF_err','SCF_jkn',...
        'Upsilon_mean','Upsilon_err','Upsilon_jkn',...
        'eta_Upsilon_mean','eta_Upsilon_err','eta_Upsilon_jkn');
end
