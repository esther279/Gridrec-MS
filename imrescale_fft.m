function img_rescale = imrescale_fft(img, scale)
    %% subpixel precision rescaling 
    if scale == 1 || isnan(scale)
        img_rescale = img;
        return
    end

    if scale > 1.2 || scale < 0.8
        warning('Too large FFT scaling can make issues ')
    end
    
    N = size(img);

    % fmat = fftshift(fft(eye(N)));
    for i = 1:2
        if N(1) ~= N(2) || i == 1
            ind_x = real(zeros(N(i),1, 'like', img));
            ind_x(:) = (0:N(i)-1)/N(i)-0.5; 
            grid = ind_x*ind_x';
            grid = -2i*pi*N(i)/scale*grid;
            W{i} = exp(grid)'/N(i);  % matrix of fourier transformation 
        else
            W{2} = W{1};
        end
    end    
    
    % fftshift can be theoretically avoided but I did not find how
    fimg = fftshift(fft2(fftshift(img)));
    
    if size(img,3) > 1
        % process large number of images together 
        fimg2 = W{1}*reshape(fimg,N,[]);
        fimg2 = reshape(fimg2, N,N,[]);
        fimg2 = permute(fimg2, [2,1,3]);
        fimg2 = reshape(fimg2, N,[]);
        img_rescale = (W{1}*fimg2);  % rescale and fft back 
        img_rescale = reshape(img_rescale, N,N,[]);
        img_rescale = permute(img_rescale, [2,1,3]);
    else
        fimg2 = W{1}*fimg;
        img_rescale = (W{2}*fimg2.').';
    end
    
    

end

