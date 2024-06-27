function [r,p,theta,omega] = mxy_snapshot_extract(snapfile,vars,system,lattice)
%MXY_SNAPSHOT_EXTRACT Extracts snapshot values r, p, te, om
%   INPUT:
%   snapfile    location of snapshot file
%   vars        defines which parameters should be extracted
%               values:     'r'     for positions
%                           'p'     for momenta
%                           't'     for spin angles
%                           'o'     for spin momenta
%               example:    'rp'    spatial data
%                           'rpto'  all data
%   system      the system to be studied
%               values:     'xy'    static xy model
%                           'mxy'   mobile xy model
%                           'vm'    Vicsek model
%                           'fvm'   frozen Vicsek model
%   lattice     lattice required for static system.
%               values:     's'     square lattice
%                           't'     trigonal lattice
%   OUTPUT:
%   r           positions in R^(2 x N)
%   p           momenta in R^(2 x N)
%   theta       angles in R^(1 x N)
%   omega       spin momenta in R^(1 x N)

    %% initialization
    r=double.empty;
    p=double.empty;
    theta=double.empty;
    omega=double.empty;
    
    %% check if anything has to be done at all
    if (contains(vars,'r') || contains(vars,'p') || ...
            contains(vars,'t') || contains(vars,'o'))
        
        snap=importdata(snapfile);
        N = length(snap(1,:));
        %% theta values
        if (contains(vars,'t'))
            theta = snap(1,:);
        end
        
        %% omega values
        if (contains(vars,'o'))
            omega = snap(2,:);
        end
        
        %% rvalues
        %  Since the line containing the r variable contains both rx and ry, 
        %  it is split into two lines by matlab.
        if (strcmp(system,"mxy") || strcmp(system,"vm"))
            if (contains(vars,'r'))
                if (strcmp(system,"mxy") )
                    i_r = 3;
                elseif ( strcmp(system,"vm"))
                    i_r = 2;
                end
                r = zeros(2,N);
                r(1,1:N/2)=snap(i_r,1:2:end);
                r(1,N/2+1:end)=snap(i_r+1,1:2:end);
                r(2,1:N/2)=snap(i_r,2:2:end);
                r(2,N/2+1:end)=snap(i_r+1,2:2:end);
            end

            %% p values
            %  Same as with r
            if (contains(vars,'p'))
                p = zeros(2,N);
                p(1,1:N/2)=snap(5,1:2:end);
                p(1,N/2+1:end)=snap(6,1:2:end);
                p(2,1:N/2)=snap(5,2:2:end);
                p(2,N/2+1:end)=snap(6,2:2:end);    
            end
        elseif (( strcmp(system,"xy") || strcmp(system,"fvm") )  ...
                && lattice == 't' && contains(vars,'r'))
            r = zeros(2,N);
            for i = 0 :sqrt(N)/2 -1
                r(1,sqrt(N)* 2 * i + 1 :sqrt(N)*(2 * i+1))=0:sqrt(N)-1;
                r(1,sqrt(N)* (2 * i + 1) + 1 :sqrt(N)*(2 * i+2))=.5 + [0:sqrt(N)-1];
                r(2,sqrt(N)* 2 * i + 1 :sqrt(N)*(2 * i+1))=sqrt(3)/2* (2 * i);
                r(2,sqrt(N)* (2 * i + 1) + 1 :sqrt(N)*(2 * i+2))=sqrt(3)/2* (2 * i + 1);
            end
        elseif (( strcmp(system,"xy") || strcmp(system,"fvm") ) ...
                && lattice == 's' && contains(vars,'r'))
            r = zeros(2,N);
            for i = 0 :sqrt(N)/2 -1
                r(1,sqrt(N)* 2 * i + 1 :sqrt(N)*(2 * i+1))=0:sqrt(N)-1;
                r(1,sqrt(N)* (2 * i + 1) + 1 :sqrt(N)*(2 * i+2))=[0:sqrt(N)-1];
                r(2,sqrt(N)* 2 * i + 1 :sqrt(N)*(2 * i+1))=(2 * i);
                r(2,sqrt(N)* (2 * i + 1) + 1 :sqrt(N)*(2 * i+2))=(2 * i + 1);
            end
        end
    else
        fprintf('Unknown mode "%s"\n', vars);
    end

end

