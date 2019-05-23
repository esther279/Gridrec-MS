function [x_crop, y_crop] = fun_crop_images(x, y)

% Crops (from the edge, ~align the center) the images to have the same size

dim_x = size(x);
cen_x = round(size(x)/2);
dim_y = size(y);
cen_y = round(size(y)/2);

dim = [min(dim_x(1), dim_y(1))  min(dim_x(2), dim_y(2))];
half = floor(dim/2);

x_crop = x(cen_x(1)-half(1)+1:cen_x(1)+half(1), cen_x(2)-half(2)+1:cen_x(2)+half(2));
y_crop = y(cen_y(1)-half(1)+1:cen_y(1)+half(1), cen_y(2)-half(2)+1:cen_y(2)+half(2));


end

