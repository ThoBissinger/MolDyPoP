function [t,y] = combine_cf(tlog,ylog,tlin,ylin)
%COMBINE_CF Combines correlation functions at logarithmic and linear times
%   simple: determines the time where logarithmic sampling becomes more
%   coarse-grained than linear sampling. Uses logarithmic data on short
%   time and linear data on long time. Maybe too simple, but not super
%   wrong
%   tlog       short-time interval with logarithmic sampling
%   ylog       function values on tlog
%   tlin       long-time interval with linear sampling
%   ylin       function values on tlin
    logseps=tlog(2:end)-tlog(1:end-1);
    linsep=tlin(2)-tlin(1);
    i_switch_log=find(logseps > linsep,1);
    if isempty(i_switch_log)
        i_switch_lin=find(tlin > max(tlog));
        if isempty(i_switch_lin)
            t = tlog;
            y = ylog;
        else
            t = [tlog,tlin(i_switch_lin:end)];
            y = [ylog,ylin(i_switch_lin:end)];
        end
    else
        t_switch=tlog(i_switch_log);
        i_switch_lin=find(tlin > t_switch, 1);
        t = [tlog(1:i_switch_log),tlin(i_switch_lin:end)];
        y = [ylog(1:i_switch_log),ylin(i_switch_lin:end)];
    end
end

