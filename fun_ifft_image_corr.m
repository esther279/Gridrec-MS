function [img_recon_0,  img_recon_corr] = fun_ifft_image_corr(cartesianGridInterpolatedFFT, interpolationCorrectionMatrix)

% Invert FFT to get the recon. image
img_recon_0 = fftshift(ifft2(ifftshift(cartesianGridInterpolatedFFT)));


% Apply interpolationCorrectionMatrix
img_recon_corr = real(img_recon_0);
interpolationCorrectionMatrix = interpolationCorrectionMatrix(1:size(img_recon_corr,1),1:size(img_recon_corr,2));

interpolationCorrectionMatrix(interpolationCorrectionMatrix>0.2*mean(interpolationCorrectionMatrix(:)))=0;
%interpolationCorrectionMatrix(interpolationCorrectionMatrix>mean(interpolationCorrectionMatrix(:)))=0;


img_recon_corr = img_recon_corr.*interpolationCorrectionMatrix;
% img_recon_corr = real(img_recon_0); only for testing

end

