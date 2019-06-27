function resolution = fun_aligned_FSC(img1, img2, param)
% resolution = fun_aligned_FSC(img1, img2, param)
% Aligns images and computes the Fourier shell correlation
% by calling the function 'fourier_shell_corr_3D_2'

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Alignment parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
align_images = param.align_images;
remove_ramp = 0;     % Try to remove ramp from whole image before initial alignment
image_prop = 'complex'; % = 'complex' or = 'phasor' (phase with unit amplitude) or = 'phase'  (Note: phase should not be used if there is phase wrapping)
cropy = [];             % Custom image cropping e.g. [100:200]- Default half size of the probe
cropx = [];             % Custom image cropping e.g. [100:200]- Default half size of the probe
flipped_images = 0; % If images are taken with a horizontal flip, e.g. 0 & 180 for tomography
guessx = [];       % Some initial guess for x alignment
guessy = [];

%%%%%%%%%%%%%%%%%%%%%%
%%% FSC parameters %%%
%%%%%%%%%%%%%%%%%%%%%%
taper = param.taper;             % Pixels to taper images - Increase until the FSC does not change anymore
SNRt = param.SNRt;          % SNRt = 0.2071 for 1/2 bit threshold for average of 2 images
                        % SNRt = 0.5 for 1 bit threshold for average of 2 images
thickring = param.thickring;

freq_thr = 0.05;        % When estimating the resolution, ignore intersections below this freq threshold
FS = 12; % font size

% img1 = recons(1).object;
% img2 = recons(2).object;

if flipped_images
    img2 = fliplr(img2);
end


% Crop images - default is half the size of the probe on each side plus


if ~isempty(guessx)
    switch sign(guessx)
        case 1
            img1 = img1(:,1+guessx:end);
            img2 = img2(:,1:end-guessx);
        case -1
            img1 = img1(:,1:end+guessx);
            img2 = img2(:,1-guessx:end);
    end
end
if ~isempty(guessy)
    switch sign(guessy)
        case 1
            img1 = img1(1+guessy:end,:);
            img2 = img2(1:end-guessy,:);
        case -1
            img1 = img1(1:end+guessy,:);
            img2 = img2(1-guessy:end,:);
    end
end

% Remove ramp 
if remove_ramp
    display('Removing ramp for initial alignment')
    img1 = remove_linearphase_v2(img1,ones(size(img1)),100);
    img2 = remove_linearphase_v2(img2,ones(size(img2)),100);
end
    

%% Alignment %%%
if align_images
    display(sprintf('\nInitial alignment'))
    switch lower(image_prop)
        case  'complex'
            imgalign1 = img1;
            imgalign2 = img2;
            display('Registering complex valued images')
        case 'phasor'
            imgalign1 = ones(size(img1)).*exp(1i*angle(img1));
            imgalign2 = ones(size(img1)).*exp(1i*angle(img2));
            display('Registering phasor of complex valued images')
        case 'phase'
            imgalign1 = angle(img1);
            imgalign2 = angle(img2);
            display('Registering phase of complex valued images')
    end

    %%% Initial alignment %%%
    upsamp = 100;
    displ = 10;
    W = 1;
    x1 = [];%[1:150];
    x2 = x1;
    y1 = [];%[1:238];
    y2 = y1;
    % imgalign2 = shiftpp2(imgalign2,10,-10); % To test range adjustment
    [subim1, subim2, delta, deltafine, regionsout] = registersubimages_2(imgalign1,imgalign2, x1, y1, x2, y2, upsamp, displ,1);
    delta_coarse = delta;
    %%% Fine alignment (second round) %%%

    % Remove ramp for fine alignment
    display('Removing ramp for fine alignment')
    %%% A patch for deltafine large
    if max(regionsout.y2+round(delta(1)))>size(img2,1)
        warning('First subpixel registration refinement found large values')
        regionsout.y2 = [min(regionsout.y2):size(img2,1)-round(delta(1))];
        regionsout.y1 = regionsout.y2;
    end
    if max(regionsout.x2+round(delta(2)))>size(img2,2)
        warning('First subpixel registration refinement found large values')
        regionsout.x2 = [min(regionsout.x2):size(img2,2)-round(delta(2))];
        regionsout.x1 = regionsout.x2;
    end
    %%%
    subimg1 = img1(regionsout.y1,regionsout.x1);
    subimg2 = img2(regionsout.y2+round(delta(1)),regionsout.x2+round(delta(2)));
    if remove_ramp
        subimg1 = remove_linearphase_v2(subimg1,ones(size(subimg1)),100);
        subimg2 = remove_linearphase_v2(subimg2,ones(size(subimg2)),100);
    end
    switch lower(image_prop)
        case  'complex'
            subimgalign1 = subimg1;
            subimgalign2 = subimg2;
            display('Registering complex valued images')
        case 'phasor'
            subimgalign1 = ones(size(subimg1)).*exp(1i*angle(subimg1));
            subimgalign2 = ones(size(subimg1)).*exp(1i*angle(subimg2));
            display('Registering phasor of complex valued images')
        case 'phase'
            subimgalign1 = angle(subimg1);
            subimgalign2 = angle(subimg2);
            display('Registering phase of complex valued images')
    end

    % Fine alignment %
    display(sprintf('\nFine alignment'))
    [subim1, subim2, delta, deltafine, regionsout] = registersubimages_2(subimgalign1,subimgalign2, x1, y1, x2, y2, upsamp, displ,1);
end

subim1 = img1 ;
subim2 = img2 ;

%% Tapering %%%
filterx = fract_hanning_pad(size(subim1,2),size(subim1,2),size(subim1,2)-2*taper);
filterx = fftshift(filterx(1,:));
filterx = repmat(filterx,[size(subim1,1) 1]);
filtery = fract_hanning_pad(size(subim1,1),size(subim1,1),size(subim1,1)-2*taper);
filtery = fftshift(filtery(:,1));
filtery = repmat(filtery,[1 size(subim1,2)]);
filterxy = filterx.*filtery;

% Taper subimages %
subim1 = subim1.*filterxy;% + (1-filterxy).*mean(subim1(:));
subim2 = subim2.*filterxy;% + (1-filterxy).*mean(subim2(:));

%% Computing the FSC
[FSC T freq] = fourier_shell_corr_3D_2(subim1,subim2,0,SNRt, thickring);

if 1
    freq_fine = 0:5e-4:max(freq);
    freq_fine_normal = freq_fine/max(freq);
    FSC_fine = interpn(freq, FSC, freq_fine, 'spline');
    T_fine = interpn(freq, T, freq_fine, 'spline');
    idx_intersect = abs(FSC_fine-T_fine)<5e-4;
    intersect_array = FSC_fine(idx_intersect);
    range = freq_fine_normal(idx_intersect);
    if length(range)<2
        range = [0 1];
        intersect_array = [1 1];
    end

    if param.plot 
        f = figure(param.plot); subplot(2,3,4);
        plot(freq_fine_normal, FSC_fine, 'linewidth',2);
        hold on; plot(freq_fine_normal, T_fine, 'r','linewidth',2);
        hold on; plot(range, intersect_array, 'go','markersize',6,'MarkerFaceColor','none','linewidth',2);
        set(gca,'fontweight','bold','fontsize',FS,'xtick',[0:0.1:1],'ytick',[0:0.1:1]);
        grid on; axis([0 1 0 1]);
        if SNRt == 0.2071 
            legend('FSC','1/2 bit');
        elseif SNRt == 0.5
            legend('FSC','1 bit');
        else
            legend('FSC',num2str(SNRt));
        end
        xlabel('Spatial frequency/Nyquist')
        %set(gcf,'Outerposition',[0 450 1000 700])    %[left, bottom, width, height
        set(gcf,'color','w');
    end

    pixel_size = param.pixelsize; %recons(1).p.dx_spec(1)*1e9; % nm
    search_range = 1; % to ignore the crossings before freq_thr
    try
        iii=2;
        while search_range
            if range(iii)>freq_thr
                range_start = range(iii); 
                search_range = 0;
            else
                iii = iii+1;
            end
        end
        resolution = [pixel_size/range_start, pixel_size/range(end)];
        catch
        resolution = [0 0];
    end
end

% res_1  = resolution(1);
% res_2 = resolution(2);



%% plot positions
% figure(11);
% plot(recons(1).p.positions(:,1),recons(1).p.positions(:,2),'ko');
% hold on
% plot(recons(2).p.positions(:,1),recons(2).p.positions(:,2),'ro');
% hold off
% 
% figure(22);
% plot(recons(1).p.positions(:,2),'k');
% hold on
% plot(recons(2).p.positions(:,2),'r');
% hold off
