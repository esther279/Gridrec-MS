function x_out = fun_padarray_2D(x, pad_dim, val)
% x_out = fun_padarray_2D(x, pad_dim, val)
% pad_dim(1,2,3,4) = (top, bottom, left, right)
%
% Example:
% x =
% 
%      1     2
%      3     4
% x_out = fun_padarray_2D(x, [2 3 1 4], [100 101 0 50])
% 
% x_out =
% 
%      0   100   100    50    50    50    50
%      0   100   100    50    50    50    50
%      0     1     2    50    50    50    50
%      0     3     4    50    50    50    50
%      0   101   101    50    50    50    50
%      0   101   101    50    50    50    50
%      0   101   101    50    50    50    50

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



x_12 = [ones([pad_dim(1), size(x,2)])*val(1); x; ones([pad_dim(2), size(x,2)])*val(2);];

x_temp = x_12';
x_1234 = [ones([pad_dim(3), size(x_temp,2)])*val(3);  x_temp; ones([pad_dim(4), size(x_temp,2)])*val(4);];

x_out = x_1234';


