
figure

n_angle=60;
angle_res = 2 * pi / n_angle;
cmapsize = n_angle;

cmap=hsv(cmapsize);
rgbcmp(1:cmapsize,1)=[linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3)];
rgbcmp(1:cmapsize,2)=[zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3),linspace(1,0.05,cmapsize/3)];
rgbcmp(1:cmapsize,3)=[linspace(1,0.05,cmapsize/3),zeros(1,cmapsize/3),linspace(0,0.95,cmapsize/3)];




outerRadius = 1;
innerRadius = 0;
numberOfSectors = 256;
graylevel = 0;
% Get polar coordinates of each point in the domain
[x, y] = meshgrid(-outerRadius : outerRadius);
[theta, rho] = cart2pol(x, y); % theta is an image here.

% Set up color wheel in hsv space.
hueImage = (theta + pi) / (2 * pi);     % Hue is in the range 0 to 1.
hueImage = ceil(hueImage * numberOfSectors) / numberOfSectors;   % Quantize hue 
saturationImage = ones(size(hueImage));      % Saturation (chroma) = 1 to be fully vivid.

% Make it have the wheel shape.
% Make a mask 1 in the wheel, and 0 outside the wheel.
wheelMaskImage = rho >= innerRadius & rho <= outerRadius;
% Hue and Saturation must be zero outside the wheel to get gray.
hueImage(~wheelMaskImage) = 0;
saturationImage(~wheelMaskImage) = 0;
% Value image must be 1 inside the wheel, and the normalized gray level outside the wheel.
normalizedGrayLevel = grayLevel / 255;
valueImage = ones(size(hueImage)); % Initialize to all 1
valueImage(~wheelMaskImage) = normalizedGrayLevel;	% Outside the wheel = the normalized gray level.

% Combine separate h, s, and v channels into a single 3D hsv image.
hsvImage = cat(3, hueImage, saturationImage, valueImage);
% Convert to rgb space for display.
rgb = hsv2rgb(hsvImage);
% Flip left to right to make it more like the usual CIE color arrangement you see.
rgb = fliplr(rgb);  % Note I think fliplr works only with color images in R2016a and later.

% Display the final color wheel.
imshow(rgb);
% Add a box around it with the pixel coordinates.
% Comment out the line below if you don't want that.
axis off;
 
