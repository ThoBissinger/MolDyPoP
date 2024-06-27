function [Phi] = NelsonFisher_Phi(y,eta)
%NELSONFISHER_PHI Eq (3.13b) in Nelson-Fisher 1977
%   Takes in y and eta and returns Phi_eta
    Phi = (y <= 1) + ...
        (y > 1) .* (y + sqrt(y.^2 - 1)).^(-eta);
end
% Example for running this script
% 
% r_vals=linspace(1e-1,148,5e2);
% t_vals=linspace(0,1e3,1e3);
% r_mat=kron(r_vals,ones(size(t_vals')));
% t_mat=kron(t_vals,ones(size(r_vals')))';
% c=1e-2;
% eta=.25;
% y = c * t_mat ./ r_mat;
% Phi = NelsonFisher_Phi(y,eta);