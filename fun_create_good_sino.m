function [sino, pad_dim] = fun_create_good_sino(img_temp, angle_singleslice, Np)

for ii = 1:length(angle_singleslice)
    theta = angle_singleslice(ii);
    img_rotate = imrotate(flipud(img_temp),theta,'bilinear','crop'); % bicubic
    sino(:,ii) = sum(img_rotate,1);
end

arbitrary_constant = Np/2;
temp = round(sqrt(sum(size(img_temp).^2))) + arbitrary_constant; 

if mod(temp,2)==0
   temp = temp+1; % Odd so that it has a center for rotation 
end

pad_dim(1) = ceil((temp - size(sino,1))/2);
pad_dim(2) = floor((temp - size(sino,1))/2);