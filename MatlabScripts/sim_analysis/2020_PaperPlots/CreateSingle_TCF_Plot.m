function [h_cf_plot,h_fit_plot] = CreateSingle_TCF_Plot(t,cf,qbin,i_N,i_T,i_q)
%CREATESINGLE_TCF_PLOT Plots TCF function and fit for parameter choice
%   Parameters
%   t           time
%   cf          correlation function (full cell)
%   i_N         index of N
%   i_T         index of T
%   i_q         index of q
%   Calling
%   Load data into matfile, e.g.
%   >> data=matfile('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/mxy/rho_3.00_dynamics_LinearTime.mat');
%   and then run e.g.
%   >> i_N=5; i_T=10; i_q = 2; t=data.averaging_times; cf=data.gmperpmperp; qbin=data.qbin; close all; [h1,h2]=CreateSingle_TCF_Plot(t,cf,qbin,i_N,i_T,i_q);
    fitfunc_DO=@(times,corrfunc_1,c) corrfunc_1 * exp(-c(1) * times/2) .* (cos(c(2) * times) + .5*c(1)/c(2)*sin(c(2) * times));
    t_vals=t{i_N,i_T};

    N_q=length(qbin{i_N,i_T});
    N_t=length(t_vals);
    q_indices = (i_q):N_q:N_t*N_q; 
    cf_vals=real(cf{i_N,i_T}(q_indices));
    
    c=fit_DampedOscillator_RealSpace(t_vals,cf_vals,10,1,'omega_1');
    h_cf_plot=plot(t_vals,cf_vals); hold on; 
    h_fit_plot=plot(t_vals,fitfunc_DO(t_vals,cf_vals(1),c));
end

