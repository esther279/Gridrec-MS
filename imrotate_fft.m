function img = imrotate_fft(img, angle, axis )
    
    Np = size(img);
    if nargin  < 3
        axis = 3;
    end
    
    if axis == 1
       img = permute(img, [3,2,1]);
    elseif axis == 2
       img = permute(img, [1,3,2]);
    end
    

    shift = 2*tan(angle/180*pi);  % Im not sure why 2, but it works as imshift from matlab
    shift_grid = linspace(-shift, shift, Np(1))' * Np(1)/4;

       
    xgrid = (fftshift((0:Np(2)-1)/Np(2))-0.5);
    if builtin( 'isa', img, 'gpuArray' )
        shift_grid = gpuArray(single(shift_grid));
        xgrid=gpuArray(single(xgrid));
    end
    X =  exp((-2i*pi)*(shift_grid(:)*xgrid));



    if angle > 0
        img = ifft(bsxfun(@times, fft(img,[],1),X'),[],1);
        img = ifft(bsxfun(@times, fft(img,[],2),X),[],2);
    else
         img = ifft(bsxfun(@times, fft(img,[],2),X),[],2);
         img = ifft(bsxfun(@times, fft(img,[],1),X'),[],1);
    end
        

    if axis == 1
       img = permute(img, [3,2,1]);
    elseif axis == 2
       img = permute(img, [1,3,2]);
    end

  
end
