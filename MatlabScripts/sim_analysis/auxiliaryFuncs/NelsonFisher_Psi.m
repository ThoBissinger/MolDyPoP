function [Psi] = NelsonFisher_Psi(y,eta,n_trunc)
%NELSONFISHER_PSI Eq (A10) in Nelson-Fisher 1977
%   Can probably be optimized, but I don't care here.
%   In modern definitions, I think that NF's hypergeometric functions of
%   the form F(a,b,c;d) would be written as F(a,b;c;d) in modern
%   definitions, so I use that here, as well.
    ind_small_y = find(y <= 1);
    ind_large_y = find(y > 1);
    small_y = y(ind_small_y);
    large_y = y(ind_large_y);
    Psi=zeros(size(y));
    if ( numel(ind_small_y) > 0 )
        for n = 0:n_trunc
            Psi(ind_small_y) = Psi(ind_small_y) ...
                + small_y.^(2 * n) / ((2*n+1)^2 - eta^2) ...
                * gamma(n + 3/2 - eta/2) / factorial(2*n+1) ...
                .* hypergeom([n + 3/2 - eta/2, (n + 3/2 - eta/2)], 2*n+2,small_y.^2);
        end
        Psi(ind_small_y) = 2^(4-eta)*cos(pi/2*eta)*eta^2 * Psi(ind_small_y);
    end
    if numel(ind_small_y) > 0
        for n = 0:n_trunc
            Psi(ind_large_y) = Psi(ind_large_y) ...
                + (-1)^(n+1) * gamma(n + 1/2 - eta/2) / gamma(n + 3/2 + eta/2) ...
                .* hypergeom([n + 3/2 - eta/2, -n + 1/2 - eta/2], 1,1./large_y.^2);
        end
        Psi(ind_large_y) = 2^(2-eta)*pi*eta^2 ./ large_y.^(3-eta).* Psi(ind_large_y);
        Psi(ind_large_y) = Psi(ind_large_y) ...
            + 2^(3-eta)*pi^2*eta^2 * sin(pi/2*eta)/(gamma(1+eta)*sin(pi*eta)) ...
            * (large_y.^2 - 1).^(eta - 1) ./ large_y.^(1 + eta);
    end
end

