function out = fract_hanning(outputdim,unmodsize);
% fract_hanning(outputdim,unmodsize)
% out = Square array containing a fractional separable Hanning window with
% DC in upper left corner.
% outputdim = size of the output array
% unmodsize = Size of the central array containing no modulation.
% Creates a square hanning window if unmodsize = 0 (or ommited), otherwise the output array 
% will contain an array of ones in the center and cosine modulation on the
% edges, the array of ones will have DC in upper left corner.

% Manuel Guizar - February 8, 2007

% Slight update on August 17, 2009
% Added a warning
% Manuel Guizar

if nargin > 2,
    error('Too many input arguments'),
elseif nargin == 1,
    unmodsize = 0;
end

if outputdim < unmodsize,
    error('Output dimension must be smaller or equal to size of unmodulated window'),
end

if unmodsize<0,
    unmodsize = 0;
    warning('Specified unmodsize<0, setting unmodsize = 0')
end

N = [0:outputdim-1];
% N = ifftshift([-floor(outputdim/2):ceil(outputdim/2)-1]);
% N = [-floor(outputdim/2):ceil(outputdim/2)-1];
[Nc,Nr] = meshgrid(N,N);

if unmodsize == 0,
    out = (1+cos(2*pi*Nc/outputdim)).*(1+cos(2*pi*Nr/outputdim))/4;
else
    % Columns modulation
    out = (1+cos(2*pi*(Nc- floor((unmodsize-1)/2) )/(outputdim+1-unmodsize)))/2;
    if floor((unmodsize-1)/2)>0,
        out(:,1:floor((unmodsize-1)/2)) = 1;
    end
    out(:,floor((unmodsize-1)/2) + outputdim+3-unmodsize:length(N)) = 1;
    % Row modulation
    out2 = (1+cos(2*pi*(Nr- floor((unmodsize-1)/2) )/(outputdim+1-unmodsize)))/2;
    if floor((unmodsize-1)/2)>0,
        out2(1:floor((unmodsize-1)/2),:) = 1;
    end
    out2(floor((unmodsize-1)/2) + outputdim+3-unmodsize:length(N),:) = 1;
    
    out = out.*out2;
end  
% out = ifftshift(out);
% one-edge at Nc = floor((unmodsize-1)/2)
% other-edge at Nc = floor((unmodsize-1)/2) + (outputdim+1-unmodsize)
%%% FINISH UP THIS CODE TO DO THE LOW PASS RECONSTRUCTION

return;
end


