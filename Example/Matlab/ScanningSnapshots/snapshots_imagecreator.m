% clear all
close all

% nsnap = 400;
nsnap = 373;
% sqrtN=32; L = 17; rho=3.55;
% sqrtN=64; L = 34; rho=3.55;
sqrtN=256; L = 136; rho=3.55;
N=sqrtN^2;

T=.13;
Tstr='.13';

runnr=7;

fold=sprintf("/data/scc/thobi/200513_MXYModel/rho_%.2f/sqrtN_%d/T_%s/run_%d/output/",rho,sqrtN,Tstr,runnr);
% fold="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/src_eclipse/output";

t=0;
dt = .01;
dt_rate = 1.025714;
% dt = 0.01; dt_rate = 1.028724; N = 200; (1 - dt_rate^N)/(1 - dt_rate)*dt

angle_res=pi/48;


read_switch = 0;
load_switch = 0;
mapcreate_switch = 0;
if (read_switch)
    rx_vals=zeros(nsnap,N);
    ry_vals=zeros(nsnap,N);
    px_vals=zeros(nsnap,N);
    py_vals=zeros(nsnap,N);
    theta_vals=zeros(nsnap,N);
    w_vals=zeros(nsnap,N);
    for i = 1:nsnap
        file=sprintf("%s/%s%d%s",fold,"snapshot_integ_",i-1,".out");
        snap=importdata(file);
        theta_vals(i,:)=snap(1,:);
        w_vals(i,:)=snap(2,:);
        rx_vals(i,1:N/2)=snap(3,1:2:end);
        rx_vals(i,N/2+1:end)=snap(4,1:2:end);
        ry_vals(i,1:N/2)=snap(3,2:2:end);
        ry_vals(i,N/2+1:end)=snap(4,2:2:end);
        px_vals(i,1:N/2)=snap(5,1:2:end);
        px_vals(i,N/2+1:end)=snap(6,1:2:end);
        py_vals(i,1:N/2)=snap(5,2:2:end);
        py_vals(i,N/2+1:end)=snap(6,2:2:end);
    end

    temperatures = (sum((px_vals.^2 + py_vals.^2 + w_vals.^2)')/3/N)';
    dirs = atan2(py_vals,px_vals);
    save('snap_data')
end
if (load_switch)
    load('snap_data')
end

% cmapsize=60;
if (mapcreate_switch)
    n_angle = 2 * pi / angle_res;
    cmapsize = n_angle;
    wrap_theta_vals=wrapToPi(theta_vals);
    sx = .75*cos(theta_vals);
    sy = .75*sin(theta_vals);
    colormap(jet);
    c_map = hsv(n_angle);
    cmap=hsv(cmapsize);
    rgbcmp(1:cmapsize,1)=[linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3)];
    rgbcmp(1:cmapsize,2)=[zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3)];
    rgbcmp(1:cmapsize,3)=[linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3)];
end
% bins=linspace(0,2*pi,cmapsize);
% angles=mod(angles,2*pi);
% col=discretize(angles,bins);

i = 372;
figure(1)
for j = 1 : n_angle
    angle = - pi + (j - 1) * angle_res;
    indices=find((angle < wrap_theta_vals(i,:)).*(wrap_theta_vals(i,:) < angle + angle_res));
    quiver(rx_vals(i,indices),ry_vals(i,indices),...
        sx(i,indices),sy(i,indices),0,...
        'color',rgbcmp(j,:));

    axis([0 L 0 L])
    hold on;
end
hold off;
%     title({sprintf('t = %.2f', dt * i), 'J = 1, U = 1'});
% title({sprintf('NVE Simulation of MXY model at t = %.2f', t), sprintf('N = %d^2, rho = %.2f, T = %.2f, U = J = 1',sqrtN,rho,T)});
  % set the seconds per image
% open the video writer
