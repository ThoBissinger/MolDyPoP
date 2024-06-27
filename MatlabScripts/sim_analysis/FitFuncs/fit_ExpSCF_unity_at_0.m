function [c] = fit_ExpSCF(rbin,SCF,r_min,r_max)
%FIT_EXPSCF Fits data (rbin,SCF) to an exponential function
%   Designed for fitting spin SCFs.
%   INPUT
%   rbin            r values
%   SCF             spin correlation function
%   r_min           starting value for r (we want long-range agreement, 
%                   short range can be ignored)
%   r_max           maximum fit distance
%   OUTPUT
%   c               fit object
    

    f_EXP = fittype('exp(-x/xi)');
    if (r_min > 0)
        index_select = find(rbin > r_min);
        rbin = rbin(index_select);
        SCF = SCF(index_select);
    end
    if (r_max > 0)
        index_select = find(rbin < r_max);
        rbin = rbin(index_select);
        SCF = SCF(index_select);
    end
%     weights=ones(size(SCF));
    weights=rbin.^2; % We want to fit large values best.
%     weights(1)=.1; % First entry not that important, possible finite size problems
%     cutoff_vec=find(SCF<threshold);
%     trust_end=length(SCF);
%     if(~ isempty(cutoff_vec))
%         trust_end=max(cutoff_vec(1)-1,1);
%         weights(trust_end+1:end) = 0.1;
%     end

%     init_corrlength = 1 / (log(SCF(end)/SCF(1)) / ( rbin(1) - rbin(end)))
    init_corrlength = 1;
%     c0=[init_factor,init_corrlength,0];
    c0=[init_corrlength];
%     fitop = fitoptions('Weights',weights,...
%         'lower',[0, 0, 0], 'upper', [1, Inf, 0]);
    c = fit(rbin(:), SCF(:), f_EXP,'StartPoint',c0,'Weights',weights,...
        'Lower',[0], 'Upper', [Inf]);
%     errfunc=@(c) norm( (fitfunc_EXP(c) - SCF) .* weights);
    
%     options = optimset('MaxFunEvals',4000);
%     c=fminsearch(errfunc,c0,options);
end

