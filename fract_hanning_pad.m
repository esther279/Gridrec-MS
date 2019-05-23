function out = fract_hanning_pad(outputdim,filterdim,unmodsize);
% fract_hanning_pad(outputdim,filterdim,unmodsize)
% out = Square array containing a fractional separable Hanning window with
% DC in upper left corner.
% outputdim = size of the output array
% filterdim = size of filter (it will zero pad if filterdim<outputdim
% unmodsize = Size of the central array containing no modulation.
% Creates a square hanning window if unmodsize = 0 (or ommited), otherwise the output array 
% will contain an array of ones in the center and cosine modulation on the
% edges, the array of ones will have DC in upper left corner.

% Manuel Guizar - August 17, 2009

if nargin > 3,
    error('Too many input arguments'),
elseif nargin == 1,
    unmodsize = 0;
    filterdim = outputdim;
end

if outputdim < unmodsize,
    error('Output dimension must be smaller or equal to size of unmodulated window'),
end

if outputdim < filterdim,
    error('Filter cannot be larger than output size'),
end

if unmodsize<0,
    unmodsize = 0;
    warning('Specified unmodsize<0, setting unmodsize = 0')
end

out = zeros(outputdim);
out(round(outputdim/2+1-filterdim/2):round(outputdim/2+1+filterdim/2-1),...
    round(outputdim/2+1-filterdim/2):round(outputdim/2+1+filterdim/2-1)) ...
    = fftshift(fract_hanning(filterdim,unmodsize));
out = fftshift(out);

return;
end

