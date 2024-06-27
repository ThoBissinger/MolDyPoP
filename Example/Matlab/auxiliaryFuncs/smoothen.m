function [y] = smoothen(x, weight,recurrence_depth)
%SMOOTHEN smoothens a function by setting the value x_i at i equal to 
%  y_i = (x_(i-1) + x_(i+1) + weight * x_i) / weight + 2.
%  Iterates recurrence_depth times
if(recurrence_depth > 1)
    y = ([0, x(1), x] + [x,x(end),0] + weight*[0,x,0]) / ( 2 + weight);
    y = y(2:end-1);
    y = smoothen(y,weight,recurrence_depth - 1);
else
    y = x;
end
end

