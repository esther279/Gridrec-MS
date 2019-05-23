function [lookupTableOfConvolventInFourierSpace, numberOfSupportNeighbours, tblspcg] = fun_interpolation(C, N_freq)

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

