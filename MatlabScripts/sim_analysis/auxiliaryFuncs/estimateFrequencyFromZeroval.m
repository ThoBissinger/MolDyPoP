function omega = estimateFrequencyFromZeroval(t,y)
%ESTIMATEFREQUENCYFROMZEROVAL Estimates frequency of an oscillating signal
%from its first zero value.
%   Uses linear interoplation to find the zero value and then computes
%   frequency as pi/2/t. Assumes the signal y to be a nonvanishing function
%   times a cosine with cos(omega*t).
    transindex=find(y < 0,1); % Finds the index of the zero transition
    if ( isempty(transindex) ) % Will not lead to sensible results, but nothing in this case will.
        transtime = t(end);
    else
        if (abs(y(transindex)) < 1e-5) % To avoid numerical problems in the 
                              % linear interpolation, we don't fit to data 
                              % too close to 0.
           transtime = t(transindex);
        elseif (abs(y(transindex-1)) < 1e-5)
           transtime = t(transindex-1);
        else % Linear interpolation
           transtime=t(transindex-1) + ...
               1/(y(transindex)/y(transindex-1)-1)...
               *(t(transindex)-t(transindex-1));
        end
    end
    omega=pi/2/transtime;
end

