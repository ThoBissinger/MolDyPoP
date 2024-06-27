function [t_new,y_new] = continue_cf(t,y,n,t_width)
%CONTINUE_CF Continues a correlation function
%   This function assumes that the correlation function just is scaled 
%   periodic in the sense that it contains a periodic part that gets damped
%   over time. That is, the whole function can be recovered by knowing one
%   period and how much the function decreases over that period.
%   Here, the continuation is carried out by identifying all periods in a
%   given function (via assuming that each period ends in a peak) and then
%   gluing the function to itself
%   t       time interval over which the function is sampled
%   y       function values on t
%   n       number of times the function should be continued
%   t_width t width for peak fit
    
    dt = (t(end) - t(1)) / ( numel(t) - 1);
    dy = (y(2:end) - y(1:end-1)) / dt;

    ind_peak = find( (dy(2:end) .* dy(1:end-1) < 0) ) + 1;
        % Peaks should be visible in signe changes. There could be multiple
        % values with sign changes around a peak, but due to the quadratic 
        % fit we do around the peak, it is enough that the actual peak is
        % within t_width of the peaks found here

    i_p = ind_peak(end); % Assuming we really just hit peaks with our procedure
    i_fit = find((t > t(i_p) - t_width) .* (t < t(i_p) + t_width)); % range for the fit
    t_fit = t(i_fit);
    y_fit = y(i_fit);

    p = polyfit(t_fit,y_fit,2);
    t_p = - p(2)/2/p(1);
    y_p = p(1)*t_p^2 + p(2)*t_p + p(3);
%     tt = linspace(t_fit(1),t_fit(end),1e4);
%     yy = 
%     if p_vals(1) > 0
% 
%     else
% 
%     end

%     y_lastpeak = y(ind_peak(end));
%     t_lastpeak = t(ind_peak(end));

    y = y(t <= t_p);
    t = t(t <= t_p);
    
    n_t = numel(t);
    t_new = zeros(1,n_t * n);
    y_new = zeros(1,n_t * n);

    scalefac = 1;
    t_start = 0;

    for i = 1:n
        t_new((i-1)*n_t + 1: i*n_t) = t + t_start;
        y_new((i-1)*n_t + 1: i*n_t) = y * scalefac;

        scalefac = scalefac * y_p / y(1);
        t_start = t_start + t_p;
    end


end

