function interpolationCorrectionMatrix = fun_interpolation_corr_matrix(C, lookupTableOfConvolventInFourierSpace, alpha)

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

