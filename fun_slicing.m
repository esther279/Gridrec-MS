function [sino_ms, img_thick] = fun_slicing(img, N_layer)

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

