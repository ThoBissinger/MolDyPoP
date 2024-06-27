function [S] = NelsonFisher_S_rt(r,c,t,eta)
%NELSONFISHER_S_RT Eq (3.13a) in Nelson-Fisher 1977
%   Takes in r, t, c and eta and returns S(r,t) as in (3.13a)
%   Meant to be used for t being fixed (then one can Fourier transform to a
%   correlation function C_m(q,t) in space)

%     Direct implementation, maybe useful to avoid infinities, but seems to
%     work anyway
%     S = (c * t < r) ./ r.^eta + ...
%         (c * t > r) .* (c * t + sqrt((c*t).^2 - r.^2)).^(-eta);
    S = r.^(-eta) .* NelsonFisher_Phi(c*t./r,eta);
end

