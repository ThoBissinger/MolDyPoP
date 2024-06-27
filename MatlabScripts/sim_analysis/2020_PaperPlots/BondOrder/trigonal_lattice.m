function [r,Lfinal] = trigonal_lattice(L,a)
%TRIGONAL_LATTICE Returns a trigonal lattice in a 2d box with spacing a
%   The box may not be perfectly square even if the size of L is. This is
%   due to the different lattice spacings in x and y direction for a
%   trigonal lattice.
%   To ensure periodicity, the lattice may also be a bit bigger in y
%   direction
%   Input values
%   L        length, can be a double or a 2d vector
%   a        lattice spacing, double
%   Output values
%   r        trigonal lattice positions, 2xN with N the number of particles
%   Lfinal   actual periodic lattice size
    if numel(L) == 1
        L = [L;L];
    elseif numel(L) > 2
        fprintf("Warning: in function trigonal_lattice:\n   More than two elements in L. Only using the first 2\n");
    end
    dx = a;
    dy = sqrt(3)/2*a;
    Nx = floor((L(1)-dx/2)/dx);
    rx = (0:Nx)*dx;
    Ny = floor(L(2)/dy);
    if (rem(Ny,2) == 1) % There has to be an even number of rows, otherwise the resulting crystal is not periodic
        Ny = Ny + 1;
    end
    ry = (0:Ny)*dy;
    N = Nx * Ny;
    rx_offsets=.5*a*rem(0:Ny,2);
    r = zeros(2,numel(rx),numel(ry));
    for i_y = 1:Ny
        r(1,:,i_y) = rx + rx_offsets(i_y);
        r(2,:,i_y) = ry(i_y);
    end
    r = reshape(r,2,[]);
    Lfinal = [max(r(1,:))+.5*a;max(r(2,:))+sqrt(3)/2*a];
end

