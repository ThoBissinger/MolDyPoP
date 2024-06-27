function [S] = StructureFactormxy(r,q)
%STRUCTUREFACTORMXY Calculates structure factor for mxy data
%   INPUT
% %   Lx      Length of sim box in x direction
% %   Ly      Length of sim box in y direction
%   r       particle positions in R^(2 x N)
%   q       wave vectors in R^(2 x m)
%   OUTPUT
%   S       Structure factor at wave vectors q in R^m

N = length(r(1,:));
m = length(q(1,:));
S = zeros(1,m);
for i = 1:m
    S(i) = abs( sum( exp(1i * (q(1,i)*r(1,:) + q(2,i)*r(2,:)) ) ) )^2;
end
S = S/N;
end

