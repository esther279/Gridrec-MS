function img = imshift_fft(img, x,y, apply_fft)
    %% subpixel precision shifting of images 
    if nargin  < 3
        y = x(:,2);
        x = x(:,1);
    end
    if nargin < 4
        apply_fft = true;
    end
    
    Np = size(img);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [nr,nc]=size(img);
%     Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
%     Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
%     [Nc,Nr] = meshgrid(Nc,Nr);
%     img2 = ifft2(fft2(img).*exp(2i*pi*(y*Nr/nr+x*Nc/nc)));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if apply_fft
         img = fft2(img);
    end

    
    % process large number of images together 
    xgrid = (fftshift((0:Np(2)-1)/Np(2))-0.5);
    X = reshape((x(:)*xgrid)',1,Np(2),[]);
    X =  exp((-2i*pi)*X);
    img = bsxfun(@times, img,X);
    ygrid = (fftshift((0:Np(1)-1)/Np(1))-0.5);
    Y = reshape((y(:)*ygrid)',Np(1),1,[]);
    Y =  exp((-2i*pi)*Y);
    img = bsxfun(@times, img,Y);
        
    if apply_fft
        img = ifft2(img);
    end
  
end
