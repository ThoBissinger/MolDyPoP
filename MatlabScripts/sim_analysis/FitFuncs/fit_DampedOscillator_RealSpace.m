function [c] = fit_DampedOscillator_RealSpace(times,corrfunc,n_period,weightexp,omega_choice)
%FIT_DAMPEDOSCILLATOR_REALSPACE Fits data (times,corrfunc) to an exponential function
%   This is meant for (time) correlation functions, so corrfunc starts with
%   a positive value. n_period is the number of periods that should be
%   included. weightexp is the exponential factor that is associated
%   with the dependence of the weight of a certain value t with its width,
%   i.e. w(i) = (t(i+1)-t(i))^weightexp. 
%   PARAMETERS:
%   times               time values
%   corrfunc            function values of the correlation function
%   n_period >= 1,      seems to work well: n_period = 4
%   weightexp >= 0,     seems to work well: weightexp \in [0,1]. 
%   omega_choice        which omega should be used. Either omega_0 or
%                       omega_1

%% Begin
    i_neg=find(corrfunc<=0,1);
    % considering the first time an index turns negative as a quarter period
    if (~isempty(i_neg))
        period_time = 4*times(i_neg);
        c0=[1,.5*pi/times(i_neg)];
    else
        period_time = times(end);
        c0=[1,0];
    end
        i_multi_period=find(times< n_period * period_time,1,'last'); 
        if(i_multi_period >= length(times))
            i_multi_period = length(times)-1; % To avoid indexing conflicts.
        end
        fitfunc_DO=@(c) corrfunc(1)*exp(-c(1) * times/2) .* (cos(c(2) * times) + .5*c(1)/c(2)*sin(c(2) * times));
        fitfunc_DO_NoTR=@(c) corrfunc(1)*exp(-c(1) * times/2) .* cos(c(2) * times) ;

        weights=ones(size(times));
%         weights(1:i_multi_period) = (times(2:i_multi_period+1) - times(1:i_multi_period)).^weightexp;
%         weights(i_multi_period+1:end) = 0;

        errfunc=@(c) norm( (fitfunc_DO(c) - corrfunc) .* weights);
        errfunc_NoTR=@(c) norm( (fitfunc_DO_NoTR(c) - corrfunc) .* weights);
%         options = optimset('MaxFunEvals',1e3*numel(times),'MaxIter',1e2*numel(times));
        options = optimset('MaxFunEvals',1e2*numel(times));
        c=fminsearch(errfunc_NoTR,c0,options); % Fit less error prone without TR term (1/om is tricky)
        c=fminsearch(errfunc,c,options);

        if (omega_choice == "omega_0")
            fitfunc_DO=@(c) corrfunc(1)*exp(-c(1) * times/2) .* ...
                (cos(sqrt(c(2)^2 - .25*c(1)^2)  * times) + ...
                .5*c(1)/sqrt(c(2)^2 - .25*c(1)^2)*sin(sqrt(c(2)^2 - .25*c(1)^2) * times)) * ...
                ( c(1) < 1e10) + ...
                corrfunc(1)*exp(-c(1) * times/2);
            errfunc=@(c) norm( (fitfunc_DO(c) - corrfunc) .* weights);
            c0 = [c(1),sqrt(c(2)^2 + .25 * c(1)^2)];
            c=fminsearch(errfunc,c0,options);
        end
end

