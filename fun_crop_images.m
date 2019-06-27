function [x_crop, y_crop] = fun_crop_images(x, y)
% [x_crop, y_crop] = fun_crop_images(x, y)
% Crops (from the edge, ~align the center) the images to have the same size

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


dim_x = size(x);
cen_x = round(size(x)/2);
dim_y = size(y);
cen_y = round(size(y)/2);

dim = [min(dim_x(1), dim_y(1))  min(dim_x(2), dim_y(2))];
half = floor(dim/2);

x_crop = x(cen_x(1)-half(1)+1:cen_x(1)+half(1), cen_x(2)-half(2)+1:cen_x(2)+half(2));
y_crop = y(cen_y(1)-half(1)+1:cen_y(1)+half(1), cen_y(2)-half(2)+1:cen_y(2)+half(2));


end

