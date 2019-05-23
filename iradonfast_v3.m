function [img,H] = iradonfast_v3(varargin)
%IRADON Inverse Radon transform.
%   I = iradon(R,THETA) reconstructs the image I from projection data in the
%   2-D array R.  The columns of R are parallel beam projection data.
%   IRADON assumes that the center of rotation is the center point of the
%   projections, which is defined as ceil(size(R,1)/2).
%
%   THETA describes the angles (in degrees) at which the projections were
%   taken.  It can be either a vector containing the angles or a scalar
%   specifying D_theta, the incremental angle between projections. If THETA
%   is a vector, it must contain angles with equal spacing between them.  If
%   THETA is a scalar specifying D_theta, the projections are taken at
%   angles THETA = m * D_theta; m = 0,1,2,...,size(R,2)-1.  If the input is
%   the empty matrix ([]), D_theta defaults to 180/size(R,2).
%
%   IRADON uses the filtered backprojection algorithm to perform the inverse
%   Radon transform.  The filter is designed directly in the frequency
%   domain and then multiplied by the FFT of the projections.  The
%   projections are zero-padded to a power of 2 before filtering to prevent
%   spatial domain aliasing and to speed up the FFT.
%
%   I = IRADON(R,THETA,INTERPOLATION,FILTER,FREQUENCY_SCALING,OUTPUT_SIZE)
%   specifies parameters to use in the inverse Radon transform.  You can
%   specify any combination of the last four arguments.  IRADON uses default
%   values for any of these arguments that you omit.
%
%   INTERPOLATION specifies the type of interpolation to use in the
%   backprojection. The default is linear interpolation. Available methods
%   are:
%
%      'nearest' - nearest neighbor interpolation 
%      'linear'  - linear interpolation (default)
%      'spline'  - spline interpolation
%      'pchip'   - shape-preserving piecewise cubic interpolation
%      'cubic'   - same as 'pchip'
%      'v5cubic' - the cubic interpolation from MATLAB 5, which does not
%                  extrapolate and uses 'spline' if X is not equally spaced.
%
%   FILTER specifies the filter to use for frequency domain filtering.
%   FILTER is a string that specifies any of the following standard filters:
% 
%   'Ram-Lak'     The cropped Ram-Lak or ramp filter (default).  The    
%                 frequency response of this filter is |f|.  Because this
%                 filter is sensitive to noise in the projections, one of
%                 the filters listed below may be preferable.
%   'Shepp-Logan' The Shepp-Logan filter multiplies the Ram-Lak filter by
%                 a sinc function.
%   'Cosine'      The cosine filter multiplies the Ram-Lak filter by a  
%                 cosine function.
%   'Hamming'     The Hamming filter multiplies the Ram-Lak filter by a
%                 Hamming window.
%   'Hann'        The Hann filter multiplies the Ram-Lak filter by a 
%                 Hann window.
%   'parzen'      The parzen filter multiplies the Ram-Lak filter by a 
%                 Parzen window.  Guizar - Nov 30 2010
%   
%   FREQUENCY_SCALING is a scalar in the range (0,1] that modifies the
%   filter by rescaling its frequency axis.  The default is 1.  If
%   FREQUENCY_SCALING is less than 1, the filter is compressed to fit into
%   the frequency range [0,FREQUENCY_SCALING], in normalized frequencies;
%   all frequencies above FREQUENCY_SCALING are set to 0.
% 
%   OUTPUT_SIZE is a scalar that specifies the number of rows and columns in
%   the reconstructed image.  If OUTPUT_SIZE is not specified, the size is
%   determined from the length of the projections:
%
%       OUTPUT_SIZE = 2*floor(size(R,1)/(2*sqrt(2)))
%
%   If you specify OUTPUT_SIZE, IRADON reconstructs a smaller or larger
%   portion of the image, but does not change the scaling of the data.
% 
%   If the projections were calculated with the RADON function, the
%   reconstructed image may not be the same size as the original image.
%
%   [I,H] = iradon(...) returns the frequency response of the filter in the
%   vector H.
%
%   Class Support
%   -------------
%   R can be double or single. All other numeric input arguments must be double.
%   I has the same class as R. H is double.
%
%   Example
%   -------
%       P = phantom(128);
%       R = radon(P,0:179);
%       I = iradon(R,0:179,'nearest','Hann');
%       figure, imshow(P), figure, imshow(I);
%
%   See also FAN2PARA, FANBEAM, IFANBEAM, PARA2FAN, PHANTOM, RADON.
%
%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2013/10/31 13:29:18 $
%
%   References: 
%      A. C. Kak, Malcolm Slaney, "Principles of Computerized Tomographic
%      Imaging", IEEE Press 1988.
%
% 'derivative'  - Optional input argument to use the filter for input
%    derivative of projections - Manuel Guizar - Nov 30 2010
%
% Weights for uneven angles
% If the angles are between 0 and 180 degrees it computes custom weights to
% the projections in order to allow for not equal angular sampling. Other
% angles are not considered because it would need significantly more
% checks.
% If you don't want to use this functionality, add 360 to your angles
% (theta+360)
% Manuel Guizar - Oct 10 2015

[p,theta,filter,d,interp,N,derivative] = parse_inputs(varargin{:});

% use Matlab or C-code for the linear interpolation
use_original_matlab_code = 0;
determine_weights = 1;

%%% Determine weights for uneven angular sampling %%%
if any(theta<0)
    warning('There are some theta < 0 angles. Using constant angular sampling code.')
    determine_weights = 0;
end
if any(theta>=pi)
    warning('There are some theta >= 180 angles. Using constant angular sampling code.')
    determine_weights = 0;
end
if any(diff(theta)==0)
    warning('There are some repeated angles. Using constant angular sampling code.')
    determine_weights = 0;
end

if abs(max(theta)-min(theta)-pi) > 5*mean(diff(sort(theta)))
    warning('Missing wedge is to large for weighting')
    determine_weights = 0;
end

if determine_weights
    [theta,ind_sort] = sort(theta);
    weights = theta*0;
    weights(2:end-1) = - theta(1:end-2)/2 + theta(3:end)/2;
    weights(1) = - (-pi + theta(end))/2 + theta(2)/2;
    weights(end) = - theta(end-1)/2 + (pi+theta(1))/2;
    p = p(:,ind_sort,:);
    p = bsxfun(@times, p, weights/2);
else
    p = p*pi/(2*length(theta));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = single(p);  % single will make it ~2x faster 

    
% Design the filter
len=size(p,1);   
H = designFilter(filter, len, d, derivative);
p(length(H),1)=0;  % Zero pad projections 

% In the code below, I continuously reuse the array p so as to
% save memory.  This makes it harder to read, but the comments
% explain what is going on.


p = fft(p);    % p holds fft of projections
p = bsxfun(@times, p, H); % frequency domain filtering

p = real(ifft(p));     % p is the filtered projections
p(len+1:end,:,:) = [];   % Truncate the filtered projections

% Define the x & y axes for the reconstructed image so that the origin     
% (center) is in the spot which RADON would choose.                        
center = floor((N + 1)/2);                                                 
xleft = -center + 1;                                                       
x = (1:N) - 1 + xleft;                                                     
x = repmat(x, N, 1);                                                       
                                                                          
ytop = center - 1;                                                         
y = (N:-1:1).' - N + ytop;                                                 
y = repmat(y, 1, N);           

if ((~strcmp(interp, 'linear')) || (use_original_matlab_code))
    costheta = cos(theta);
    sintheta = sin(theta);
    img = zeros(N,class(p));        % Allocate memory for the image.
end

ctrIdx = ceil(len/2);     % index of the center of the projections

% Zero pad the projections to size 1+2*ceil(N/sqrt(2)) if this
% quantity is greater than the length of the projections
imgDiag = 2*ceil(N/sqrt(2))+1;  % largest distance through image.
Nlayers = size(p,3); 

if size(p,1) < imgDiag 
   rz = imgDiag - size(p,1);  % how many rows of zeros
   pad = zeros([ceil(rz/2),size(p,2),Nlayers], 'single'); 
   p = [pad; p; pad];
   ctrIdx = ctrIdx+ceil(rz/2);
end

% Backprojection - vectorized in (x,y), looping over theta
switch interp
    case 'nearest neighbor'

      for i=1:length(theta)   
          proj = p(:,i);
          t = round(x*costheta(i) + y*sintheta(i));
          img = img + proj(t+ctrIdx);
      end
   
  case 'linear'
    if (use_original_matlab_code)
        for i=1:length(theta)  
            proj = p(:,i);
            t = x.*costheta(i) + y.*sintheta(i); 
            a = floor(t);  
            img = img + (t-a).*proj(a+1+ctrIdx) + (a+1-t).*proj(a+ctrIdx);
    %         imagesc(img);drawnow;
        end
    else
        try
            gpuDevice;
            if Nlayers == 1; warning('Using only single slice is ineffecient !!! '); end
            img = iradon_cuda( p, theta, x, y );
        catch
            warning('GPU failed, running slower CPU version')
            keyboard
            img = zeros(N,N,Nlayers);
            p = double(p);
            for i = 1:size(p,3)
                progressbar(i,Nlayers)
                img(:,:,i) = iradon_c( p(:,:,i), theta, x, y );
            end
            img = single(img);
        end
    end
    
    
  case {'spline','pchip','cubic','v5cubic'}

    interp_method = sprintf('*%s',interp); % Add asterisk to assert
                                           % even-spacing of taxis
    
    for i=1:length(theta)
        proj = p(:,i);
        taxis = (1:size(p,1)) - ctrIdx;
        t = x.*costheta(i) + y.*sintheta(i);
        projContrib = interp1(taxis,proj,t(:),interp_method);
        img = img + reshape(projContrib,N,N);
    end
    
end

% img = img*pi/(2*length(theta));


%%%
%%%  Sub-Function:  designFilter
%%%

function filt = designFilter(filter, len, d, derivative)
% Returns the Fourier Transform of the filter which will be 
% used to filter the projections
%
% INPUT ARGS:   filter - either the string specifying the filter 
%               len    - the length of the projections
%               d      - the fraction of frequencies below the nyquist
%                        which we want to pass
%
% OUTPUT ARGS:  filt   - the filter to use on the projections


order = max(64,2^nextpow2(2*len));

% First create a ramp filter - go up to the next highest
% power of 2.
if derivative
    filt = 0*( 0:(order/2) )+1;
else
    filt = 2*( 0:(order/2) )./order;
end
w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist 

switch filter
case 'ram-lak'
   % Do nothing
case 'shepp-logan'
   % be careful not to divide by 0:
   filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
case 'cosine'
   filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
case 'hamming'  
   filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
case 'hann'
   filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
case 'parzen'
   aux = parzenwin(round(2*size(filt,2)*d)-1)';
   aux = aux(round(size(aux,2)/2):round(size(aux,2)));
   filt(1:size(aux,2)) = filt(1:size(aux,2)).*aux;
   filt(size(aux,2)+1:end) = 0;
otherwise
   eid = sprintf('Images:%s:invalidFilter',mfilename);
   msg = 'Invalid filter selected.';
   error(eid,'%s',msg);
end

filt(w>pi*d) = 0;                      % Crop the frequency response
if derivative
    filt = [filt' ; -filt(end-1:-1:2)']/(1i*pi);    % Symmetry of the filter
else
    filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter
end


%%%
%%%  Sub-Function:  parse_inputs
%%%

function [p,theta,filter,d,interp,N,derivative] = parse_inputs(varargin)
%  Parse the input arguments and retun things
%
%  Inputs:   varargin -   Cell array containing all of the actual inputs
%
%  Outputs:  p        -   Projection data
%            theta    -   the angles at which the projections were taken
%            filter   -   string specifying filter or the actual filter
%            d        -   a scalar specifying normalized freq. at which to crop 
%                         the frequency response of the filter
%            interp   -   the type of interpolation to use
%            N        -   The size of the reconstructed image

if nargin<2
   eid = sprintf('Images:%s:tooFewInputs',mfilename);
   msg = 'Invalid input arguments.';
   error(eid,'%s',msg);
end

p = varargin{1};
theta = pi*varargin{2}/180;

% Default values
N = 0;                 % Size of the reconstructed image
d = 1;                 % Defaults to no cropping of filters frequency response
filter = 'ram-lak';    % The ramp filter is the default
interp = 'linear';     % default interpolation is linear
string_args = {'nearest neighbor', 'linear', 'spline', 'pchip', 'cubic', 'v5cubic', ...
      'ram-lak','shepp-logan','cosine','hamming', 'hann','parzen','derivative'};

for i=3:nargin
   arg = varargin{i};
   if ischar(arg)
      idx = strmatch(lower(arg),string_args);
      if isempty(idx)
         eid = sprintf('Images:%s:unknownInputString',mfilename);
         msg = sprintf('Unknown input string: %s.', arg);
         error(eid,'%s',msg);
      elseif numel(idx) > 1
         eid = sprintf('Images:%s:ambiguousInputString',mfilename);
         msg = sprintf('Ambiguous input string: %s.', arg);
         error(eid,'%s',msg);
      elseif numel(idx) == 1
         if idx <= 6   % It is the interpolation
            interp = string_args{idx};
         elseif (idx > 6) && (idx <= 12)
            filter = string_args{idx};
         elseif idx == 13
            derivative = true; % Input is a derivative of sinogram
         end
      end
   elseif numel(arg)==1
      if arg <=1
         d = arg;
      else
         N = arg;
      end
   else
      eid = sprintf('Images:%s:invalidInputParameters',mfilename);
      msg = 'Invalid input parameters';
      error(eid,'%s',msg);
   end
end

% If the user didn't specify the size of the reconstruction, so 
% deduce it from the length of projections
if N==0    
   N = 2*floor( size(p,1)/(2*sqrt(2)) );  % This doesn't always jive with RADON
end

% for empty theta, choose an intelligent default delta-theta
if isempty(theta)
   theta = pi / size(p,2);
end

% If the user passed in delta-theta, build the vector of theta values
if numel(theta)==1
   theta = (0:(size(p,2)-1))* theta;
end

if length(theta) ~= size(p,2)
   eid = sprintf('Images:%s:thetaNotMatchingProjectionNumber',mfilename);
   msg = 'THETA does not match the number of projections.';
   error(eid,'%s',msg);
end

if ~exist('derivative')
    derivative = false;
end

function img = iradon_cuda( sinogram_full, theta, x, y )
    %% preprocess data for the ASTRA toolbox wrapper

    [Wsin, Nproj, Nlayers] = size(sinogram_full); 
    [Nx, Ny] = size(x);

    %% create data and geometry
    Nangles = length(theta);
    assert(Nproj == Nangles, 'Wrong input size')

    vectors_all = zeros(Nangles, 12);
    for i = 1:Nangles
      % ray direction
      vectors_all(i,1) = sin(theta(i));
      vectors_all(i,2) = -cos(theta(i));
      vectors_all(i,3) = 0;	
      vectors_all(i,1:3) =  vectors_all(i,1:3);
      % center of detector
      vectors_all(i,4:6) = 0;
      % vector from detector pixel (0,0) to (0,1)  
      vectors_all(i,7) = cos(theta(i));
      vectors_all(i,8) = sin(theta(i));
      vectors_all(i,9) = 0;
      % vector from detector pixel (0,0) to (1,0) 
      vectors_all(i,10) = 0;
      vectors_all(i,11) = 0;
      vectors_all(i,12) = 1;     
    end
    %% astra settings 
    cfg.iVolX = Nx; 
    cfg.iVolY = Ny; 
    cfg.iVolZ = Nlayers; 
    cfg.iProjAngles = Nangles; 
    cfg.iProjU = Wsin; 
    cfg.iProjV = Nlayers;
    cfg.iRaysPerDet = 1; 
    cfg.iRaysPerDetDim = 1; 
    cfg.iRaysPerVoxelDim = 1; 

 


    %% call ASTRA wrapper 
    gpu = gpuDevice; 
    AvailableMemory = gpu.AvailableMemory; 
    % required memory is 2x dataset size, max data size allowed by GPU is 1024MB
    % max number of GPU projections is 1024 (limit in ASTRA code, can be fixed ... )
    Nangular_slices =  ceil(Nangles/1024); % limitation in the ASTRA code
    Nslices = ceil(numel(sinogram_full)*8 / min(1024e6, AvailableMemory)); 
    for j = 1:Nangular_slices  % split in the angular space 
        if Nangular_slices > 1 || Nslices > 1
            fprintf('Dataset does not fit to GPU memory => autosplitting\nFree memory: %iMB\tDataset size: %iMB\n\n', ceil(AvailableMemory/1e6), ceil(numel(sinogram_full)*4/1e6))
            img = zeros([Nx, Ny, Nlayers], 'single'); 
        end
        
        ind_angle = 1+(j-1)*ceil(Nangles/Nangular_slices):min(Nlayers, j*ceil(Nangles/Nangular_slices));         
        vectors = vectors_all(ind_angle,:);
        cfg.iProjAngles = length(ind_angle);
        % split the sinogram (avoid copying if possible)
        if Nangular_slices > 1
            sinogram = sinogram_full(:,ind_angle,:);
        else
            sinogram = sinogram_full; 
        end
        if Nslices == 1
            %%%% apply geometry correction to shift reconstruction into center %%%%%%%%%%%%%%%%%%%%%%%%%%
            vectors(:,4:6) = vectors(:,4:6) -(vectors(:,10:12)*cfg.iProjV/2+vectors(:,7:9)*cfg.iProjU/2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % reshape the inputs for ASTRA 
            sinogram =  reshape(sinogram,[cfg.iProjU,cfg.iProjV,cfg.iProjAngles]); 
            sinogram =  gpuArray(single(sinogram)); 
            %% ASTRA WRAPPER 
            vol = ASTRA_GPU_wrapper('bp', sinogram, cfg, vectors);
            vol = gather(vol);  % backprojection 
            if Nangular_slices > 1
                img = img + vol;
            else
                img = vol;
            end
        else  % low memory case (a bit slower) => split to several chunks rotation axis 
            for i = 1:Nslices
                progressbar(i,Nslices+1)
                ind = 1+(i-1)*ceil(Nlayers/Nslices):min(Nlayers, i*ceil(Nlayers/Nslices)); 
                p_tmp = sinogram(:,:,ind);
                cfg.iVolZ = length(ind); 
                cfg.iProjV = length(ind); 
                p_tmp =  reshape(p_tmp,[cfg.iProjU,cfg.iProjV,cfg.iProjAngles]); 
                p_tmp = gpuArray(single(p_tmp)); 
                vectors_all_tmp = vectors; 
                                
                %%%% apply geometry correction to shift reconstruction into center
                vectors_all_tmp(:,4:6) = vectors_all_tmp(:,4:6) -(vectors_all_tmp(:,10:12)*cfg.iProjV/2+vectors_all_tmp(:,7:9)*cfg.iProjU/2);
                %% ASTRA WRAPPER
                vol = ASTRA_GPU_wrapper('bp', p_tmp, cfg, vectors_all_tmp);
                vol = gather(vol);  % gather from GPU takes 30% of time !!! 
                if Nangular_slices > 1
                    img(:,:,ind) = img(:,:,ind) + vol;  % adding up to the array is also pretty slow
                else
                    img(:,:,ind) = vol;  
                end
            end
            progressbar(i,Nslices)
        end
    end
    % rotate as the output for Matlab 
    img = rot90(img, 1); 
    
    
function progressbar(n,N,w)
    % progressbar - display a progress bar
    %
    %    progressbar(n,N,w);
    %
    % displays the progress of n out of N.
    % n should start at 1.
    % w is the width of the bar (default w=20).
    %
    %   Copyright (c) Gabriel Peyr 2006

    if nargin<3
        w = 20;
    end

    % progress char
    cprog = '.';
    cprog1 = '*';
    % begining char
    cbeg = '[';
    % ending char
    cend = ']';

    p = min( floor(n/N*(w+1)), w);

    global pprev;
    if isempty(pprev)
        pprev = -1;
    end

    if not(p==pprev)
        ps = repmat(cprog, [1 w]);
        ps(1:p) = cprog1;
        ps = [cbeg ps cend];
        if n>1
            % clear previous string
            fprintf( repmat('\b', [1 length(ps)]) );
        end
        fprintf(ps);
    end
    pprev = p;
    if n==N
        fprintf('\n');
    end



