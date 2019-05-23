function [img_temp, img_orig, fimg_orig] = fun_create_good_phantom(img, zp_len)

% Padded with zero on all sides, make sure dimension is odd



img_orig = fun_padarray_2D(img, [zp_len zp_len zp_len zp_len], [0 0 0 0]);
fimg_orig = fftshift(fft2(ifftshift(img_orig)));
% if 0 
% 	recon = fftshift(ifft2(ifftshift(fimg_orig))); 
% 	figure(100); 
%     subplot(1,2,1); imagesc(log10(abs(fimg_orig))); axis equal xy tight; colorbar; title('log10(abs(fimg_orig))');
%     subplot(1,2,2); imagesc(abs(recon)); axis equal xy tight; colorbar;
% end

if mod(size(img_orig,1),2)==0 % pad the image so that the image has a center pixel
    img_temp = fun_padarray_2D(img_orig, [1 0 1 0], [0 0 0 0]);
else
    img_temp = img_orig;
end



end

