function [FSC T r n] = fourier_shell_corr_3D_2(img1,img2,dispfsc,SNRt,thickring) 
% [FSC T r] = fourier_shell_corr_3D(img1,img2,dispfsc,SNRt,thickring) 
% Computes the Fourier shell correlation between img1 and img2. It can also
% compute the threshold function T. Images can be complex-valued.
% Can handle non-cube arrays but assumes the voxel is isotropic
% dispfsc = 1;  % Display results
% SNRt          % Power SNR for threshold, popular options:
                % SNRt = 0.5; 1 bit threshold for average
                % SNRt = 0.2071; 1/2 bit threshold for average
% thickring	% Normally the pixels get assigned to the closest integer pixel ring in Fourier domain. With thickring the thickness of 
		% the rings is increased by thickring, so each ring gets more pixels and more statistics
%M. van Heela, and M. Schatzb, "Fourier shell correlation threshold
%criteria," Journal of Structural Biology 151, 250-262 (2005)
% Manuel Guizar 19-Apr-2011
% Added option of thick frequency rings. Thickness is the last parameter in
% pixels

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


% A = 3;
% img1 = rand(100,100,100);
% img2 = img1 + A*rand(100,100,100);
% img1 = img1 + A*rand(100,100,100);
% dispfsc = 1;
% SNRt = 1/A^2;

if any(size(img1) ~= size(img2))
    error('Images must be the same size')
end

if ~exist('dispfsc')
    dispfsc = 0;
end

if ~exist('thickring')
    thickring = 0;
end
[ny nx nz] = size(img1);
nmax = max([nx ny nz]);

x = ifftshift([-fix(nx/2):ceil(nx/2)-1])*floor(nmax/2)/floor(nx/2);
y = ifftshift([-fix(ny/2):ceil(ny/2)-1])*floor(nmax/2)/floor(ny/2);
if nz ~= 1
    z = ifftshift([-fix(nz/2):ceil(nz/2)-1])*floor(nmax/2)/floor(nz/2);
else
    z = 0;
end
[X Y Z] = meshgrid(single(x),single(y),single(z));
index = round(sqrt(X.^2+Y.^2+Z.^2));
clear X Y Z
rnyquist = floor(nmax/2);
r = [0:rnyquist];
%
% for ii = 40;%1:size(index,3),
% figure(1); 
% imagesc(fftshift(index(:,:,ii))),
% colorbar; axis xy
% end
%


F1 = fftn(img1);
F2 = fftn(img2);

for ii = r;
    if thickring == 0,
        auxF1 = F1(index==ii);
        auxF2 = F2(index==ii);
    else
        auxF1 = F1( (index>=(ii-thickring/2))&(index<=(ii+thickring/2)) );
        auxF2 = F2( (index>=(ii-thickring/2))&(index<=(ii+thickring/2)) );
    end
    C(ii+1)  = sum(auxF1.*conj(auxF2));
    C1(ii+1) = sum(auxF1.*conj(auxF1));
    C2(ii+1) = sum(auxF2.*conj(auxF2));
    n(ii+1) = size(auxF1,1);  % Number of points
end
%
FSC = abs(C)./(sqrt(C1.*C2));

if exist('SNRt')
    T = (  SNRt + 2*sqrt(SNRt)./sqrt(n+eps) + 1./sqrt(n)  )./...
        (  SNRt + 2*sqrt(SNRt)./sqrt(n+eps) + 1  );
end

if dispfsc
    figure(1);
    plot(r/rnyquist,FSC,'Linewidth',1);
    if exist('SNRt')
        hold on,
        plot(r/rnyquist,T,'--r','Linewidth',1);
        hold off,
        if SNRt == 0.2071,
            legend('FSC','1/2 bit threshold')
        elseif SNRt == 0.5,
            legend('FSC','1 bit threshold')
        else
            legend('FSC',['Threshold SNR = ' num2str(SNRt)]);
        end
    else
        legend('FSC')
    end
    xlim([0 1])
    ylim([0 1.1])
    xlabel('Spatial frequency/Nyquist')
    
end
