function [c] = fit_GaussOscillator(times,corrfunc)
%FIT_GAUSSOSCILLATOR Fits data (times,corrfunc) to an oscillating Gaussian
%   This is meant for (time) correlation functions, so f starts with a
%   positive value. The exponential 
    negative_indices=find(corrfunc<=0);
    i_neg=negative_indices(1);
    % considering the first time an index turns negative as a quarter period
    period=find(times<4*times(i_neg)); 
    i_period=period(end);
    fitfunc_GAUSS=@(c) corrfunc(1)*exp(-.5 * c(1) * times.^2) .* cos(c(2) * times);
    weights=zeros(size(corrfunc));
    weights(1:i_neg) = 1;
    weights(i_neg+1:i_period) = 1;
    weights(i_period+1:end) = 0;
    errfunc=@(c) norm( (fitfunc_GAUSS(c) - corrfunc) .* weights);
    c0=[1,.5*pi/negative_indices(1)];
    options = optimset('MaxFunEvals',2000);
    c=fminsearch(errfunc,c0,options);
end

