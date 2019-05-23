function [ out ] = fun_rotate( x,theta )

% https://ch.mathworks.com/matlabcentral/answers/215259-why-is-imrotate-not-rotating-an-image-about-its-center-point

Rdefault = imref2d(size(x));  
tX = mean(Rdefault.XWorldLimits);
tY = mean(Rdefault.YWorldLimits);
tTranslationToCenterAtOrigin = [1 0 0; 0 1 0; -tX -tY,1];
% theta = 45;
tRotation = [cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1];
tTranslationBackToOriginalCenter = [1 0 0; 0 1 0; tX tY,1];
tformCenteredRotation = tTranslationToCenterAtOrigin*tRotation*tTranslationBackToOriginalCenter;
tformCenteredRotation = affine2d(tformCenteredRotation);
[out,Rout] = imwarp(x,tformCenteredRotation);
%imshow(out,Rout)

end

