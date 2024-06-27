function [c] = fit_PowSCF(rbin,SCF,r_min,r_max,offset_switch)
%FIT_POWSCF Fits data (rbin,SCF) to a power function
%   Designed for fitting spin SCFs.
%   INPUT
%   rbin            r values
%   SCF             spin correlation function
%   r_min           starting value for r (we want long-range agreement, 
%                   short range can be ignored)
%   r_max           maximum fit distance
%   offset_switch   if this is 1, then there is an offset parameter in the
%                   fit function, allowing for decay to a plateau
%   OUTPUT
%   c               fit object
    if (offset_switch)
        f_POW = fittype('a*(x)^(-eta) - c');
    else
        f_POW = fittype('a*(x)^(-eta)');
%         f_POW = fittype('-eta*log(x)+log(a)');
    end
    
%     weights(1)=.1; % First entry not that important, possible finite size problems
%     weights=rbin; % For Power Law, Long-Time-Behavior is important.
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
    weights=ones(size(SCF));
%     cutoff_vec=find(SCF<threshold);
%     trust_end=length(SCF);
%     if(cutoff_vec)
%         trust_end=max(cutoff_vec(1)-1,1);
%         weights(trust_end+1:end) = 0.1;
%     end
%     init_exponent = 1 / (log(SCF(end)/SCF(1)) / ( log(rbin(1) / rbin(end))));
%     init_factor = SCF(2) / (rbin(2))^(init_exponent);
    init_exponent = .25;
    init_factor = 1;
    if (offset_switch)
        c0=[init_factor,init_exponent,0];
        Lower_lim=[0, 0, 0];
        Upper_lim=[1, Inf,1];
    else
        c0=[init_factor,init_exponent];
        Lower_lim=[0, 0.002]; % To avoid constant solution
        Upper_lim=[1, Inf];
    end
    c = fit(rbin(:), SCF(:), f_POW,'StartPoint',c0,'Weights',weights,...
        'Lower',Lower_lim, 'Upper', Upper_lim);
%     errfunc=@(c) norm( (fitfunc_EXP(c) - SCF) .* weights);
    
%     options = optimset('MaxFunEvals',4000);
%     c=fminsearch(errfunc,c0,options);
end

