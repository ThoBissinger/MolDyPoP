function [c] = fit_eta_Susceptibility_FS(chi,sqrtN_vals)
%FIT_ETA_SUSCEPTBILITY_FS Obtains eta values by a fit to the FS
%   Susceptibility chi = C sqrt(N)^(2-eta)
%   Factor C for correctioneta_vals_absM(i_T))

    chi = chi(:) ./ sqrtN_vals(:).^2;
    p = polyfit(log(sqrtN_vals(:)), log(chi(:)),1);
    c.eta=-p(1)/2;
    c.a=exp(p(2));
%  Older Implementation, using fits. Doesn't work as well as the polyfit of
%  the log
%     f_POW = fittype('a*x^(-2*eta)');
%     f_LOG = fittype('a-eta*x');
    
%     c0 = [5e-3,0];
%     Lower_lim=[0,0];
%     Upper_lim=[1,2];
%     c = fit(sqrtN_vals(:), chi(:), f_POW,'StartPoint',c0,...
%         'Lower',Lower_lim, 'Upper', Upper_lim);

%     c = fit(log(sqrtN_vals(:)), log(chi(:)), f_POW,'StartPoint',c0,...
%         'Lower',Lower_lim, 'Upper', Upper_lim);
end

