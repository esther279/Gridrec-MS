function [sino_fft, sino_ms] = fun_generate_sino_ms(angle_array, N_layer, img_temp, pad_dim)

N_theta = length(angle_array);

fprintf('----- \n');
for thetaIndex = 1:N_theta
    theta = angle_array(thetaIndex);
    img_rotate = imrotate(flipud(img_temp),theta,'bilinear','crop'); 
    img_rotate = imnoise(img_rotate, 'Gaussian',0.01);
    
    [temp_sino, thick_img] = fun_slicing(img_rotate, N_layer); 
    use = temp_sino;
    if sum(pad_dim)~=0
        sino_ms(:, :, thetaIndex) = fun_padarray_2D(use, [0 0 pad_dim], [0 0 0 0]);
    else
        sino_ms(:, :, thetaIndex) = use;
    end
    
    fimg_ms = fftshift( fft2(ifftshift(sino_ms(:, :, thetaIndex) )) );    
    fimg_ms = squeeze(fimg_ms);
    for n = 1:N_layer
        sino_fft(n, :, thetaIndex) = fimg_ms(n,:);   
    end
end


end

