function [error] = CalculateCollapseError(curves,rbins,L_vals,c,d,plotswitch)
%CALCULATEFITERROR Summary of this function goes here
%   Detailed explanation goes here
    error = 0;
    N_L = numel(L_vals);
%     rr=rbyLmin:rbyLmax;
    for i = 2 : N_L
        for j = 1: i -1
            
            rmax=min(rbins{i}(end)*L_vals(j)^c/L_vals(i)^c,rbins{j}(end));
            rmin=max(rbins{j}(1),rbins{i}(1));
            
            rr = linspace(rmin,rmax);
            error=error + mean(abs(curves{i}(rr*L_vals(i)^c/L_vals(j)^c) * L_vals(i)^d - curves{j}(rr) * L_vals(j)^d));
            % Logarithmic error. Maybe better for fitting chimq?
%             error=error + mean(abs(log(curves{i}(rr*L_vals(i)^c/L_vals(j)^c) * L_vals(i)^d) - log(curves{j}(rr) * L_vals(j)^d)));

            if (plotswitch)
                figure(10*i + 100*j)
                plot(rr,curves{i}(rr*L_vals(i)^c/L_vals(j)^c) * L_vals(i)^d,...
                    'DisplayName',sprintf('$L = (%.1f)^2$', L_vals(i)));
                hold on;
                plot(rr,curves{j}(rr) * L_vals(j)^d,...
                    'DisplayName',sprintf('$L = (%.1f)^2$', L_vals(j)));
%                 legend show
                legend('interpreter','latex');
            end
% 
            
            
            
        end
    end
end

