function [ft_vals,om_vals] = FT_correlation(t_vals,func_vals, n_points,res_func)
%FT_CORRELATION Fourier transform for correlation function
%   Calculates the Fourier transform of the data set (t_vals,func_vals)
%   using a cubic spline.
%   t_vals      time values, "x-axis"
%   func_vals   (correlation) function values, "y-axis"
%   n_points    number of data points the spline is used for
    if ( n_points <= length(t_vals) ) 
        spline_t_vals = t_vals;
        spline_func_vals = func_vals;
        n_points=length(t_vals);
    else
        spline_t_vals = linspace(t_vals(1),t_vals(end),n_points);
        spline_func_vals = spline(t_vals,func_vals,spline_t_vals);
    end
%     dt = spline_t_vals(2) - spline_t_vals(1);
    dt = (spline_t_vals(end) - spline_t_vals(1))/(n_points - 1);
    dom = 2 * pi / dt;
    
    ft_vals = 2*dt*fftshift(fft(spline_func_vals)); % possibly missing multiplication with factor 2*pi or sth
    om_vals = (-n_points/2:n_points/2-1) * dom / n_points;
    
end

