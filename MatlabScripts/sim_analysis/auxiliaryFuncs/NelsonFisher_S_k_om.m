function [S_k_om,k,om] = NelsonFisher_S_k_om(r,t,c,eta)
%NELSONFISHER_S_K_OM k-omega-space FT of Eq (3.13a) in Nelson-Fisher 1977
%   Takes in r, t, c and eta and returns S(k,om), the Fourier transform 
%   of (3.13a) in both reciprocal lattice (k) and frequency domain
%   Input
%   r       positions, has to be linearly spaced. Assumes row vector form,
%           size(r) = [1 n_r]
%   t       times, has to be linearly spaced. Assumes column vector form, 
%           size(t) = [n_t 1]
%   c       speed of sound
%   eta     critical exponent eta
%   Output
%   S       function S(k,t), Fourier transform of S(r,t) in (3.13a)
%   k       reciprocal lattice vectors
%   om      frequencies
    if (size(r,1) ~= 1)
        r = r';
    end
    if (size(t,2) ~= 1)
        t = t';
    end
    r_mat = kron(r,ones(size(t)));
    t_mat = kron(t,ones(size(r)));
    y = c*t_mat./r_mat;
    
    n_points_r=numel(r);
    dr = r(2) - r(1);
    dk = 2 * pi / dr;
    k = (-n_points_r/2:n_points_r/2-1) * dk / n_points_r;
    
    n_points_t=numel(t);
    dt = t(2) - t(1);
    dom = 2 * pi / dt;
    om = (-n_points_t/2:n_points_t/2-1) * dk / n_points_t;
    
    S_r_t = r_mat.^(-eta) .* NelsonFisher_Phi(y,eta);
    S_k_t = 2*dr*fftshift(fft(S_r_t'),1)';
    S_k_om = 2*dt*fftshift(fft(S_k_t),2);
         % FFT for every column of S_r', i.e. every row of S_r.
         % Attention: one must specify the dimension along which fftshift 
         % has to applied. We want to use the r-dimension, which is the
         % column dimension in S_r'.    
end

