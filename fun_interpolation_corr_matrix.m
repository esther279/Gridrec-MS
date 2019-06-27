function interpolationCorrectionMatrix = fun_interpolation_corr_matrix(C, lookupTableOfConvolventInFourierSpace, alpha)
% interpolationCorrectionMatrix = fun_interpolation_corr_matrix(C, lookupTableOfConvolventInFourierSpace, alpha)

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


N_freq = length(lookupTableOfConvolventInFourierSpace);
%[lookupTableOfConvolventInFourierSpace, ~, ~] = fun_interpolation(C, N_freq);
numberOfSamplingPoints = N_freq/2; % no round, but this is not an integer

scaleRatio = double(numberOfSamplingPoints)/double(numberOfSamplingPoints + 0.5);
%    Norm is only important for deriving absolute absorbtion coefficient, this is however not possible with this form of gridrec
norm = sqrt(pi/2/C/alpha); 
lookupTableOfConvolventInRealSpace = zeros(1,N_freq); %zeros(1,2*numberOfSamplingPoints);

center = ceil((N_freq+1)/2);
lookupTableOfConvolventInRealSpace(center) = norm/lookupTableOfConvolventInFourierSpace(1);
for samplingPoint = 1:floor(N_freq/2) % was 0:numberOfSamplingPoints
    if center + samplingPoint <= N_freq
        lookupTableOfConvolventInRealSpace(center + samplingPoint) = norm/lookupTableOfConvolventInFourierSpace(round(samplingPoint*scaleRatio));
    end
    lookupTableOfConvolventInRealSpace(center - samplingPoint) = norm/lookupTableOfConvolventInFourierSpace(round(samplingPoint*scaleRatio));
    % lookupTableOfConvolventInRealSpace(numberOfSamplingPoints + samplingPoint) = norm/lookupTableOfConvolventInFourierSpace(round(samplingPoint*scaleRatio)+1);
    % lookupTableOfConvolventInRealSpace(numberOfSamplingPoints+1 - (samplingPoint)) = norm/lookupTableOfConvolventInFourierSpace(round(samplingPoint*scaleRatio)+1);
end
interpolationCorrectionMatrix = lookupTableOfConvolventInRealSpace'*lookupTableOfConvolventInRealSpace;

end

