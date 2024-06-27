function [bcmean,bias,stderr,jkn] = jackknife_estimate(data,fnhand)
%JACKKNIFE_ESTIMATE Calculates the jackknife estimates for bias corrected mean and stderr
%   Given a dataset `data`, uses the typical methods for calculation of the
%   mean, the bias and the standard error via jackknife
%   fnhand is the function handle that goes into jackknife
    if size(data,1) == 1 % Data should be a column vector or a matrix
        data = data';
    end
    if nargin < 2
        fnhand = @mean;
    end
    datasize=size(data);
    n=size(data,1);

    meanval=fnhand(data);
    
    jkn=jackknife(fnhand,data);
    jkn=reshape(jkn,size(data));
    jkn_mean=mean(jkn);
    bias=(n-1)*(jkn_mean - meanval);
    bcmean=meanval -  bias;

    pseudovals=n*meanval - (n-1)*jkn;
    jackvar=1/(n*(n-1))*sum((pseudovals-bcmean).^2);
    stderr=sqrt(jackvar);
    
end

