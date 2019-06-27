function [sino_ms, img_thick] = fun_slicing(img, N_layer)
% [sino_ms, img_thick] = fun_slicing(img, N_layer)

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


N = size(img,1);
thickness = floor(N/N_layer);
thickest = thickness + mod(N,N_layer); % BAD IDEA!!
thick_edge = thickness + mod(N,N_layer)/2;

img_thick = zeros(size(img));

if N_layer>1 && 0 
    gaussFilter = gausswin(thickness*2, 6); % default ALPHA is 2.5 (the larger the narrower FWHM)
    gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.
    %gaussFilter = ones(size(gaussFilter));
    plot(gaussFilter) 

    for ii = 1:size(img,2)
        img_conv(:,ii) = conv(img(:,ii), gaussFilter,'same');
    end

    idx = 1;
    for n = 1:N_layer
        idx_array = [idx:idx+thickness-1];
        if n==ceil(N_layer/2)
            idx_array = [idx:idx+thickest-1];
        else
            idx_array = [idx:idx+thickness-1];
        end
        idx = idx_array(end)+1;
        mid = round(mean(idx_array));

        %slice = mean(img(idx_array,:),1)';
        slice = img_conv(mid,:);
        sino_ms(n,:) = slice;

        for ii = idx_array
            img_thick(ii,:) = slice';
        end
    end  
else

    idx = 1;
    for n = 1:N_layer
        if n==1 || n==N_layer
            idx_array = [idx:idx+thick_edge-1];
        else
            idx_array = [idx:idx+thickness-1];
        end
        idx = idx_array(end)+1;

        slice = sum(img(idx_array,:),1)'; % was mean
        slice = slice./thickness;
        sino_ms(n,:) = slice;

        for ii = idx_array
            img_thick(ii,:) = slice';
        end

    end  
end

end

