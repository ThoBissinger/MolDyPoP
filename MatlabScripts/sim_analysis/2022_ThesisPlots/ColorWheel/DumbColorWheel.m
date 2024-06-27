%%
cd /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/ColorWheel
figure

n_angle=1200;
angle_res = 2 * pi / n_angle;
cmapsize = n_angle;
angle_vals=-pi:angle_res:pi-angle_res/2;

cmap=hsv(cmapsize);
rgbcmp(1:cmapsize,1)=[linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3)];
rgbcmp(1:cmapsize,2)=[zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3)];
rgbcmp(1:cmapsize,3)=[linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3)];


% x=linspace(-1,1);
% [X,Y] = meshgrid(x,x);
% angles=atan2(Y,X);
% for j = 1 : n_angle
%     angle = - pi + (j - 1) * angle_res;
%     indices=find((angle < angles).*(angles < angle + angle_res));
%     plot(X(indices),Y(indices),...
%         'color',rgbcmp(j,:));
%     hold on;
% end
% 
% figure
% mesh(angles,'color',rgbcmp);
%% 
n_pix=1001;
n_pix_half = 500;
% [X,Y] = meshgrid(1:n_pix,1:n_pix);
image=zeros(n_pix,n_pix,3); %initialize
alpha=zeros(n_pix,n_pix); %initialize
% angles=atan2(Y,X);
for x = -n_pix_half:n_pix_half
    i_x = x + n_pix_half + 1;
    for y = -round(sqrt(n_pix_half^2-x^2)):1:round(sqrt(n_pix_half^2-x^2))
        cur_angle = wrapToPi(atan2(y,x) - pi/2); % messed up sth, dunno what.  Pi/2 offset corrects that
        i_y = y + n_pix_half + 1;
        [minval,c_ind]=min(abs(cur_angle-angle_vals));
        image(i_x,i_y,:)=rgbcmp(c_ind,:);
        alpha(i_x,i_y,:)=1;
    end
end

figure, imshow(image)
imwrite(image,"colorwheel.png","alpha",alpha);