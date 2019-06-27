function [Q, x_crop, y_crop] = fun_calc_error(x, y, param)
% [Q, x_crop, y_crop] = fun_calc_error(x, y, param)
% Calculate the reconstruction quality Q: correlation, rmse, FSC resolution

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


[x_crop, y_crop] = fun_crop_images(x, y);

Q(2) = sqrt(mean((x_crop(:)-y_crop(:)).^2));

temp1 = x_crop-mean(x_crop(:)); 
temp2 = y_crop-mean(y_crop(:));
Q(1) = temp1(:)'*temp2(:) / sqrt(sum(temp1(:).^2)*sum(temp2(:).^2));

Q(3:4) = fun_aligned_FSC(x_crop, y_crop, param);

% figure;
% subplot(1,2,1);
% imagesc(x_crop);
% 
% subplot(1,2,2);
% imagesc(y_crop);
