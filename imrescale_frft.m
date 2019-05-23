function img = imrescale_frft(img, scale)
    %% imrescale_fft based on fractional fourier transformation 

    %% 2d version is faster only for many stucked pictures 
    if size(img,3) == 1 || length(scale) > 1
        img = fftshift(ifft(fftshift(FRFT_Centered(img,scale(1)))));
        img = fftshift(ifft(fftshift(FRFT_Centered(img',scale(end)))))';
    else
        img = fftshift(ifft2(fftshift(FRFT_2D(img,scale(1)))));
    end

end

function Y=FRFT_Centered(X,alpha)
    % rewritten and speed / memory optimized by  Michal Odstrcil 2015, Southampton University, 
    % See A. Averbuch, "Fast and Accurate Polar Fourier Transform"

    %% it maybe works as magification lens Claus, D., & Rodenburg, J. M. (2015). Pixel size adjustment in coherent diffractive imaging within the Rayleigh–Sommerfeld regime
           
    N = size(X,1);

    grid = fftshift(-N:N-1)';

    preFactor = exp(1i*pi*grid*alpha);  % perform shift 
    Factor=exp(-1i*pi*grid.^2/N * alpha);  % propagation / scaling
    
    x_tilde=[X; zeros(size(X), class(X))];  % add oversampling 
    x_tilde= bsxfun(@times, x_tilde, Factor .* preFactor);

    XX=fft(x_tilde);
        
    YY=fft(conj(Factor));  
    
    XX = bsxfun(@times, XX,YY);
    Y=ifft(XX);
    
    Y=bsxfun(@times, Y,Factor.*preFactor);
    Y=Y(1:N,:,:);
    

end

function Y=FRFT_2D(X,alpha)
          
    N = size(X,1);
    grid = fftshift(-N:N-1)';
    [Xg,Yg] = meshgrid(grid, grid);
    preFactor = exp(1i*pi*(Xg+Yg)*alpha);  % perform shift after FFT 
    Factor=exp(-1i*pi*(Xg.^2+Yg.^2)/N(1) * alpha);  % propagation / scaling
    x_tilde = zeros(2*N, 2*N, size(X,3),  class(X));
    x_tilde(1:N, 1:N,:) = X;
    x_tilde= bsxfun(@times, x_tilde, Factor .* preFactor);
        
    XX=fft2(x_tilde);

    YY=fft2(conj(Factor));  

    XX = bsxfun(@times, XX,YY);

    Y=ifft2( XX );
    Y=bsxfun(@times, Y,Factor .* preFactor );
    Y=Y(1:N,1:N,:);

end

