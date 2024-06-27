function [fitob] = fit_DampedOscillator_FourierSpace(omega,h,h_init)
%FIT_DAMPEDOSCILLATOR_FOURIERSPACE Fits the damped oscillator in reciprocal
%frequency space. 
%   Fit contains omega_0 and gamma, and the function to be fitted is
%   2 * h_0 * gamma * omega_0^2
%        / ( (omega^2 - omega_0^2)^2 + (omega gamma)^2 )
%        + h_0 * noiselevel
%   omega_0 is the resonance frequency, gamma the damping, h_0 is the
%   value of h(t = 0), noiselevel quantifies the amount of noise in the
%   data (it is at least the value of the time spacing dt
    
    if (isrow(h))
        h = h';
    end
    if (isrow(omega))
        omega = omega';
    end
    h = h/h_init;

    [hmax, i_omega_max] = max(h);
    initguess_omega_0 = abs(omega(i_omega_max)); % h should be symmetric, but it could be that for numerical reasons, a negative omega_0 is found to be maximal
    weights = ones(size(h));
    if (initguess_omega_0 > 1e-5) 
%         initguess_gamma = 1 / (initguess_omega_0^2 * (hmax - initguess_dt));
        initguess_gamma = 2 / (hmax);
        weights(max(i_omega_max -5,1) : min(i_omega_max + 5,numel(weights))) = 1e8; % Stronger weights on the peak area.
    else
        initguess_gamma = max(omega(find(h > hmax/2)));
        initguess_omega_0 = initguess_gamma;
%         initguess_omega_0=1e-4;
%         initguess_gamma = 0.025; % Different guess for central peak
    end


%     f_DO_om = fittype('2*h_0*gamma*omega_0^2/((x^2 - omega_0^2)^2 + (x * gamma)^2) + noiselevel');
%     c0=[initguess_gamma, initguess_h_0, initguess_noiselevel, initguess_omega_0];
    f_DO_om = fittype('2*gamma*omega_0^2/((x^2 - omega_0^2)^2 + (x * gamma)^2)');
    c0=[initguess_gamma, initguess_omega_0];
%     if nargin < 3
%         LowLim=[0,0, 0, 0];
%         UpperLim=[Inf,Inf,noisemax, Inf];
%     else
%         LowLim=[0,.9*initguess_h_0, 0, 0];
%         UpperLim=[Inf,1.1*initguess_h_0,noisemax, Inf];
%     end
%     if nargin < 3
%         LowLim=[0,0, 0];
%         UpperLim=[Inf,Inf, Inf];
%     else
%         LowLim=[0,0];
%         UpperLim=[Inf,Inf];
%     end
    LowLim=[0,0];
    UpperLim=[Inf,Inf];
    fitob = fit(omega, h, f_DO_om,'StartPoint',c0,'Weights',weights,...
         'Lower',LowLim, 'Upper', UpperLim); 
end

