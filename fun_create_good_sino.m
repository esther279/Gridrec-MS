function [sino, pad_dim] = fun_create_good_sino(img_temp, angle_singleslice, Np)
% [sino, pad_dim] = fun_create_good_sino(img_temp, angle_singleslice, Np)

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



for ii = 1:length(angle_singleslice)
    theta = angle_singleslice(ii);
    img_rotate = imrotate(flipud(img_temp),theta,'bilinear','crop'); % bicubic
    sino(:,ii) = sum(img_rotate,1);
end

arbitrary_constant = Np/2;
temp = round(sqrt(sum(size(img_temp).^2))) + arbitrary_constant; 

if mod(temp,2)==0
   temp = temp+1; % Odd so that it has a center for rotation 
end

pad_dim(1) = ceil((temp - size(sino,1))/2);
pad_dim(2) = floor((temp - size(sino,1))/2);
