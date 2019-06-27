function [lookupTableOfConvolventInFourierSpace, numberOfSupportNeighbours, tblspcg] = fun_interpolation(C, N_freq)
% [lookupTableOfConvolventInFourierSpace, numberOfSupportNeighbours, tblspcg] = fun_interpolation(C, N_freq)

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


%% paramters
support_len = 2*C/pi;
half_support_len = support_len/2;
tblspcg = (N_freq/2)/half_support_len;
numberOfSupportNeighbours = ceil(half_support_len);

%% === setupInterpolationLookupTable
numberOfSamplingPoints = N_freq/2; % no round, but this is not an integer
coeff = [0.5767616E+02, 0.0, -0.8931343E+02, 0.0, ...
    0.4167596E+02, 0.0, -0.1053599E+02, 0.0, 0.1662374E+01, 0.0, ...
    -0.1780527E-00, 0.0, 0.1372983E-01, 0.0, -0.7963169E-03, 0.0, ...
    0.3593372E-04, 0.0, -0.1295941E-05, 0.0, 0.3817796E-07];
temp1 = fun_legendre(coeff, [0:numberOfSamplingPoints+1] / double(numberOfSamplingPoints) );
temp2 = fun_legendre(coeff, 0.0);
temp3 = temp1 / temp2;
%lookupTableOfConvolventInFourierSpace = padarray(temp3, [0 numberOfSamplingPoints], 0); 
lookupTableOfConvolventInFourierSpace = [temp3 zeros(1,round(numberOfSamplingPoints)+1)];  % no round, but this is not an integer


end

