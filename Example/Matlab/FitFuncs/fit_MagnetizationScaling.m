function [T_C,T_star,T_KT,crossover_M, fitfuncs] = fit_MagnetizationScaling(T,absM,sqrtN,coeff)
%FIT_MAGNETIZATIONSCALING Finite size scaling fits for the magnetization,
%   After Bramwell and Holdsworth
    if nargin < 4
        coeff=sqrt(2);
    end
    N_N = numel(absM(:,1));
    T_C = zeros(numel(N_N),1);
    T_star = zeros(numel(N_N),1);
    T_KT = zeros(numel(N_N),1);
    crossover_M = zeros(numel(N_N),1);
    fitfuncs = cell(numel(N_N),1);
    
    beta_exp=3*pi^2/128;
    
    for i_N = 1 : N_N
        crossover_M(i_N) = (1 / coeff / sqrtN(i_N))^.125;
        curM_vals=absM(i_N,:);
        i_Tstar=find(curM_vals > crossover_M(i_N),1,'last');
        T_star(i_N) = T(i_Tstar) + ...
            (T(i_Tstar+1)-T(i_Tstar)) * ...
            (crossover_M(i_N) - curM_vals(i_Tstar))/ ...
            (curM_vals(i_Tstar+1) - curM_vals(i_Tstar));
   

        fit_indices=find(T < 1.1 * T_star(i_N) ...
            & T > .5 * T_star(i_N));
        T_C0=1.5*T_star(i_N);
        
        curfitfunc=@(x,Tvar) 1/(x-T_star(i_N))^(beta_exp) * crossover_M(i_N) ...
            * (x - Tvar).^(beta_exp);%.*(Tvar<x);
        fitfuncs{i_N} = curfitfunc;
        errfunc=@(x) norm((curfitfunc(x,T(fit_indices)) - curM_vals(fit_indices))...
            .*(T(fit_indices)<x)+1e1 * (x < T_star(i_N)));
        options = optimset('MaxFunEvals',5000);
        T_C(i_N)=fminsearch(errfunc,T_C0,options);
        T_KT(i_N) = (4 * T_star(i_N) - T_C(i_N))/3;
    end
end

