function [sino_fft, sino_ms] = fun_generate_sino_ms(angle_array, N_layer, img_temp, pad_dim)
% [sino_fft, sino_ms] = fun_generate_sino_ms(angle_array, N_layer, img_temp, pad_dim)

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a
%   publication or if it is fully or partially rewritten for another
%   computing language the authors and institution should be acknowledged
%   in written form in the publication: “Data processing was carried out
%   using the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.”
%   Variations on the latter text can be incorporated upon discussion with
%   the CXS group if needed to more specifically reflect the use of the package
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%  
% This code and subroutines are part of a continuous development, they
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its
%    proper use and the correctness of the results.


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

