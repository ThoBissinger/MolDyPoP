function [q_vals,chi,chimpar,chimperp,chi_mfree,chimpar_mfree,S] = Susceptibility(r,te,L,q_max,m)
%STRUCTUREFACTORMXY Calculates structure factor for mxy data
%   INPUT
% %   Lx      Length of sim box in x direction
% %   Ly      Length of sim box in y direction
%   r       particle positions in R^(2 x N)
%   te      theta values in R^(N)
%   q       wave vectors in R^(2 x m)
%   m       spontaneous magentization    
%   OUTPUT
%   chi     Susceptibility at wave vectors q in R^m
N = length(r(1,:));
s_x = cos(te);
s_y = sin(te);
if nargin < 5
    m = 1/N*sum([s_x;s_y]');
end
e_m = m /norm(m);
te_m = atan2(m(2),m(1));
q_vals = 0:2*pi/L:q_max;
n_q = numel(q_vals);

mqx = zeros(2,n_q); % sum(s_x.* exp(1i * (q(1,i)*r(1,:))
mqy = mqx;
rhoq = mqx;
% chi = zeros(1,n_q);
for i = 1:n_q
    q=q_vals(i);
    for j = 1:2
        mqx(j,i) = sum(s_x.* exp(-1i * q*r(j,:)));
        mqy(j,i) = sum(s_y.* exp(-1i * q*r(j,:)));
        rhoq(j,i) = sum(exp(-1i * q*r(j,:)));
    end
end
mqx = mqx / sqrt(N);
mqy = mqy / sqrt(N);
rhoq = rhoq / sqrt(N);

chi = (abs(mqx(1,:)).^2 + abs(mqx(2,:)).^2 + abs(mqy(1,:)).^2 + abs(mqy(2,:)).^2)/2;
chimpar = (abs(e_m(1)*mqx(1,:)+e_m(2)*mqy(1,:)).^2 + abs(e_m(1)*mqx(2,:)+e_m(2)*mqy(2,:)).^2)/2;
chimperp = (abs(-e_m(2)*mqx(1,:)+e_m(1)*mqy(1,:)).^2 + abs(-e_m(2)*mqx(2,:)+e_m(1)*mqy(2,:)).^2)/2;

chi_mfree = (abs(mqx(1,:)-m(1)*rhoq(1,:)).^2 + abs(mqx(2,:)-m(1)*rhoq(2,:)).^2 ...
    + abs(mqy(1,:)-m(2)*rhoq(1,:)).^2 + abs(mqy(2,:)-m(2)*rhoq(2,:)).^2)/2;
chimpar_mfree = (abs(e_m(1)*mqx(1,:)+e_m(2)*mqy(1,:)-norm(m)*rhoq(1,:)).^2 ...
    + abs(e_m(1)*mqx(2,:)+e_m(2)*mqy(2,:)-norm(m)*rhoq(2,:)).^2)/2;

S = ( abs(rhoq(1,:)).^2 + abs(rhoq(2,:)).^2)/2;
% chi(i) = abs( sum( exp(1i * (q(1,i)*r(1,:) + q(2,i)*r(2,:)) ) ) )^2;
% S = S/N;
end

