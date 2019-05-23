function img = imshift_fft_ax(img, shift, ax)

    Np = size(img);
    Np(ax) = 1;
    Ng = ones(1,3);
    N = length(shift);
    Ng(ax) = N;

    grid = (fftshift((0:N-1)/N)-0.5);

    
    if ~isvector(shift)
        X = bsxfun(@times, reshape(shift,Np), reshape(grid,Ng));
        X =  exp((2i*pi)*X);
    else
        X = bsxfun(@times, shift, reshape(grid,Ng));
        X =  exp((2i*pi)*X);
    end

    img = fft(img,[],ax);


    img = bsxfun(@times, img,X);
   
    img = ifft(img,[],ax);

end

