%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Script to align images and compute FSC %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References relevant to this code:
% For using this FSC code with ptychography: J. Vila-Comamala, et al., "Characterization of high-resolution diffractive X-ray optics by ptychographic coherent diffractive imaging," Opt. Express 19, 21333-21344 (2011).
% For subpixel alignment: M. Guizar-Sicairos, et al., "Efficient subpixel image registration algorithms," Opt. Lett. 33, 156 (2008).
% For matching of phase ramp by approximate least squared error:  M. Guizar-Sicairos, et al., "Phase tomography from x-ray coherent diffractive imaging projections," Opt. Express 19, 21345-21357 (2011).
clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruction files %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openwithGUI = 0;    % Use GUI for choosing files, otherwise specify parameters below
scan1 = [31];            % Scan number of first image
scan2 = [32];            % Scan number of second image
sample_name = 'porglass_';    % File prefix
suffix = '_opt.mat';         % File suffix
analysis_folder = '../ptycho_test_data/analysis/';
% online_ptycho_S00031_S00032_600x600_test_1_c
% S00032_S00031_600x600_N2_cp_z25um_MLit300_recons
% S00031_S00032_600x600_ms_N2_z20um_test_msroi_1_3ML_it236
% S00031_S00032_600x600_ms_N3_z20um_test_1_3ML_it500
% S00031_S00032_600x600_ms_N3_z20um_test_2_3ML_it502
% S00031_S00032_600x600_ms_N2_z20um_test_3_3ML_it486
% S00031_S00032_600x600_ms_N2_z20um_test_3_3DM_it1400
% S00031_S00032_600x600_ms_N3_z20um_1_3DM_it1000
filenamewithpath1 = ['../ptycho_test_data/analysis/S00031/S00031_S00032_600x600_ms_N3_z20um_1b_3ML_it200.mat'];       % Give the full filename and path - Overrides the parameters above
filenamewithpath2 = ['../ptycho_test_data/analysis/S00032/S00031_S00032_600x600_ms_N3_z20um_1b_3ML_it200.mat'];       % Give the full filename and path - Overrides the parameters above

% x = load(file(1).filename)
% x.recon_time(end)/3600
% [x.ms_error(2,1), x.ms_error(2,end)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Alignment parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

remove_ramp = 1;     % Try to remove ramp from whole image before initial alignment
image_prop = 'phasor'; % = 'complex' or = 'phasor' (phase with unit amplitude) or = 'phase'  (Note: phase should not be used if there is phase wrapping)
cropy = [];             % Custom image cropping e.g. [100:200]- Default half size of the probe
cropx = [2000:3000];             % Custom image cropping e.g. [100:200]- Default half size of the probe
flipped_images = 0; % If images are taken with a horizontal flip, e.g. 0 & 180 for tomography
GUIguess = 0;       % To click for an initial alignment guess, ignores the values below
guessx = [];       % Some initial guess for x alignment
guessy = [];

flag_alignment = 0;
total_delta = [0.1 0.33];
% total_delta = ([0.1 0.3]+[-0.04 1.79]+[0.05 1.41])/3;  % 600x600_test_1_c
% total_delta = ([0.1 0.38]+[-0.03 1.68]+[0.05 1.24])/3; % N2_cp_z25um_MLit300_recons

%%%%%%%%%%%%%%%%%%%%%%
%%% FSC parameters %%%
%%%%%%%%%%%%%%%%%%%%%%
taper = 20;             % Pixels to taper images - Increase until the FSC does not change anymore
dispfsc = 1;            % Display FSC plot
SNRt = 0.2071;          % SNRt = 0.2071 for 1/2 bit threshold for average of 2 images
                        % SNRt = 0.5 for 1 bit threshold for average of 2 images
thickring = 2;
freq_thr = 0.05;        % When estimating the resolution, ignore intersections below this freq threshold

FS = 12; % font size

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Do not modify below %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

scanfolder1 = sprintf('S%05d', scan1(1)); % Looks for the file in this folder, I leave a variable so that the folder can be overriden
scanfolder2 = sprintf('S%05d', scan2(1));
%%% Opening file %%%
if openwithGUI
    display('Using GUI open mode')
    if exist([analysis_folder scanfolder1],'dir')
        uipath1 = [analysis_folder scanfolder1];
    else
        uipath1 = [];
    end
    if exist([analysis_folder scanfolder2],'dir')
        uipath2 = [analysis_folder scanfolder2];
    else
        uipath2 = [];
    end
    [filename, pathname] = uigetfile('*.mat','Open first reconstruction', uipath1);
    file(1).filename = fullfile(pathname,filename);
    
    [filename, pathname] = uigetfile('*.mat','Open second reconstruction', uipath2);
    file(2).filename = fullfile(pathname,filename);
else
    % Opening recons 1 %
    if ~isempty(filenamewithpath1)
        file(1).filename = filenamewithpath1;
    else
        file(1).filename = fullfile(analysis_folder,scanfolder1,[sample_name '*' suffix]);
        D = dir(file(1).filename);
        if numel(D) == 0
            error(['I did not find any file: ' file(1).filename])
        elseif numel(D) > 1
            warning(['I found many files with the mask: ' file(1).filename]);
            warning(['I selected ' D(1).name]);
        end
        file(1).filename = fullfile(analysis_folder,scanfolder1,D(1).name);
    end
    
    % Opening recons 2 %
    if ~isempty(filenamewithpath2)
        file(2).filename = filenamewithpath2;
    else
        file(2).filename = fullfile(analysis_folder,scanfolder2,[sample_name '*' suffix]);
        D = dir(file(2).filename);
        if numel(D) == 0
            error(['I did not find any file: ' file(2).filename])
        elseif numel(D) > 1
            warning(['I found many files with the mask: ' file(2).filename]);
            warning(['I selected ' D(1).name]);
        end
        file(2).filename = fullfile(analysis_folder,scanfolder2,D(1).name);    
    end
end

if exist(file(1).filename,'file')
    display(['Loading: ' file(1).filename])
    recons(1) = load(file(1).filename,'object','p');
    file(1).titlestring = strrep(file(1).filename,'_','\_');
else
    error(['Not found: ' file(1).filename])
end
if exist(file(2).filename,'file')
    display(['Loading: ' file(2).filename])
    recons(2) = load(file(2).filename,'object','p');
    file(2).titlestring = strrep(file(2).filename,'_','\_');
else
    error(['Not found: ' file(2).filename])
end

img1 = recons(1).object;
img2 = recons(2).object;

if flipped_images
    img2 = fliplr(img2);
end

asize = recons(1).p.asize;

% Show phase images (not cropped)%
figure(2)
set(gcf,'Outerposition',[1 600 500 500],'color','w')    %[left, bottom, width, height
imagesc(angle(img1));
axis xy equal tight
colormap bone
aux = angle(img1(round(asize(1)/2):end-round(asize(1)/2),round(asize(2)/2):end-round(asize(2)/2)));
caxis([min(aux(:)) max(aux(:))])
title(file(1).titlestring)
figure(3)
imagesc(angle(img2));
axis xy equal tight
colormap bone
aux = angle(img2(round(asize(1)/2):end-round(asize(1)/2),round(asize(2)/2):end-round(asize(2)/2)));
caxis([min(aux(:)) max(aux(:))])
title(file(2).titlestring)
set(gcf,'Outerposition',[500 600 500 500],'color','w')    %[left, bottom, width, height


% Crop images - default is half the size of the probe on each side plus
% whatever needed to make them of equal size
minsize1 = min(size(img1,1),size(img2,1));
minsize2 = min(size(img1,2),size(img2,2));
if isempty(cropy)
    cropy = [round(asize(1)/2):minsize1-round(asize(1)/2)];
end
if isempty(cropx)
    cropx = [round(asize(2)/2):minsize2-round(asize(1)/2)];
end
img1 = img1(cropy,cropx);
img2 = img2(cropy,cropx);

if GUIguess
    figure(2)
    display(['Click on a feature on figure 2'])
    [xin yin] = ginput(1);
    figure(3)
    display(['Click on a feature on figure 3'])
    [xin2 yin2] = ginput(1);
    guessx = round(xin-xin2);
    guessy = round(yin-yin2);
end

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
    

figure(4)
set(gcf,'Outerposition',[1 1 500 476],'color','w')    %[left, bottom, width, height
imagesc(angle(img1));
axis xy equal tight
colormap bone
title(file(1).titlestring)
figure(5)
imagesc(angle(img2));
axis xy equal tight
colormap bone
title(file(2).titlestring)
set(gcf,'Outerposition',[500 1 500 476],'color','w')    %[left, bottom, width, height

%
try
    recons(1).ms_error(1)
    recons(1).ms_error(end)
catch
end

%% Alignment %%%
if flag_alignment
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
    subimg1 = remove_linearphase_v2(subimg1,ones(size(subimg1)),100);
    subimg2 = remove_linearphase_v2(subimg2,ones(size(subimg2)),100);

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
else
    % Use existing alignement
    display(sprintf('\nUse existing alignement'))
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

    upsamp = 100;
    %errors = dftregistration(fft2(imgalign1), fft2(imgalign2), upsamp)
    temp2 = shiftpp2(imgalign2, -total_delta(1), -total_delta(2));
    %errors = dftregistration(fft2(imgalign1), fft2(temp2), upsamp)

    subim1 = imgalign1;
    subim2 = temp2;
end
errors1 = dftregistration(fft2(imgalign1), fft2(imgalign2), upsamp);
errors2 = dftregistration(fft2(subim1), fft2(subim2), upsamp);
total_delta = errors1(3:4) - errors2(3:4)

%% Trying to stitch images together, but the phase wrap are different..
if 0
    subim1 = remove_linearphase_v2(subim1,ones(size(subim1)),100);
    subim2 = remove_linearphase_v2(subim2,ones(size(subim2)),100);
    subim3 = subim2;
    subim3(:,[1:1000 2000:2973]) = subim1(:,[1:1000 2000:2973]);  
    figure(100); imagesc(angle(subim3)); colormap bone; axis xy equal tight;
end

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

figure(4)
set(gcf,'Outerposition',[1 1 500 476])    %[left, bottom, width, height
if strcmpi(image_prop,'phase')
    imagesc(subim1);
else
    imagesc(angle(subim1));
end
axis xy equal tight
colormap bone
title(file(1).titlestring)
figure(5)
if strcmpi(image_prop,'phase')
    imagesc(real(subim2));
else
    imagesc(angle(subim2));
end
axis xy equal tight
colormap bone
title(file(2).titlestring)
set(gcf,'Outerposition',[500 1 500 476])    %[left, bottom, width, height

%% Computing the FSC
[FSC T freq] = fourier_shell_corr_3D_2(subim1,subim2,0,SNRt, thickring);

display(sprintf('\nPixel size = %.3f nm',recons(1).p.dx_spec(1)*1e9))

if dispfsc
    freq_fine = 0:1e-3:max(freq);
    freq_fine_normal = freq_fine/max(freq);
    FSC_fine = interpn(freq, FSC, freq_fine, 'spline');
    T_fine = interpn(freq, T, freq_fine, 'spline');
    idx_intersect = abs(FSC_fine-T_fine)<1e-4;
    intersect_array = FSC_fine(idx_intersect);
    range = freq_fine_normal(idx_intersect);
    if length(range)<2
        range = [0 1];
        intersect_array = [1 1];
    end

    f = figure; 
    plot(freq_fine_normal, FSC_fine, 'b.-','linewidth',1);
    hold on; plot(freq_fine_normal, T_fine, 'r','linewidth',1);
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
    set(gcf,'Outerposition',[0 450 1000 700])    %[left, bottom, width, height
    set(gcf,'color','w');

    pixel_size = recons(1).p.dx_spec(1)*1e9; % nm
    search_range = 1; % to ignore the crossings before freq_thr
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
    st_title = sprintf('%s\n%s \n flipped_images %d, remove_ramp %d \n cropx [%d:%d], align %d, delta [%.3f, %.3f] \nPixel size %.2f nm\n FSC: taper %d, thickring %d, intersect (%.3f, %.3f) \n Resolution (%.2f, %.2f) nm', ...
        file(1).filename, file(2).filename, flipped_images, remove_ramp, cropx(1), cropx(end),flag_alignment,total_delta(1), total_delta(2),pixel_size, taper, thickring, range_start, range(end), pixel_size/range_start, pixel_size/range(end))
    title(st_title,'interpreter','none')
end

res_1  = resolution(1);
res_2 = resolution(2);



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
