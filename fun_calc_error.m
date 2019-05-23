function [Q, x_crop, y_crop] = fun_calc_error(x, y, param)

% Calculate the error, actually the quality Q: correlation, rmse, FSC resolution

% Crops (from the edge) the images to have the same size
% Calculates the quality (Q, 1 for true) and the RMSE between the images

% dim_x = size(x);
% cen_x = round(size(x)/2);
% dim_y = size(y);
% cen_y = round(size(y)/2);
% 
% dim = [min(dim_x(1), dim_y(1))  min(dim_x(2), dim_y(2))];
% half = round(dim/2);
% 
% x_crop = x(cen_x(1)-half(1)+1:cen_x+half(1), cen_x(2)-half(2)+1:cen_x+half(2));
% y_crop = y(cen_y(1)-half(1)+1:cen_y+half(1), cen_y(2)-half(2)+1:cen_y+half(2));

[x_crop, y_crop] = fun_crop_images(x, y);

Q(2) = sqrt(mean((x_crop(:)-y_crop(:)).^2));

temp1 = x_crop-mean(x_crop(:)); 
temp2 = y_crop-mean(y_crop(:));
Q(1) = temp1(:)'*temp2(:) / sqrt(sum(temp1(:).^2)*sum(temp2(:).^2));

Q(3:4) = fun_aligned_FSC(x_crop, y_crop, param);

% figure;
% subplot(1,2,1);
% imagesc(x_crop);
% 
% subplot(1,2,2);
% imagesc(y_crop);