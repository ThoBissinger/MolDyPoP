function [c] = fit_eta_Magnetization_FS(M,sqrtN_vals)
%FIT_ETA_MAGNETIZATION_FS Obtains eta values by a fit to the FS
%   magnetization M = C(1/sqrt(2N))^eta/2
%   Factor C for correction
%   M magnetization
    f_POW = fittype('c*(sqrt(2)*x)^(-eta/2)');
    c0 = [1,.25];
    Lower_lim=[0,0];
    Upper_lim=[Inf,4];
    c = fit(sqrtN_vals(:), M(:), f_POW,'StartPoint',c0,...
        'Lower',Lower_lim, 'Upper', Upper_lim);
end
