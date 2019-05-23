%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2017-01-19
%
% Original: see gridRec.py from TOMCAT and 
% F. Marone and M. Stampanoni, "Regridding reconstruction algorithm for 
% real-time tomographic imaging," J.Synchrotron Radiat. 19, 1029-1037 (2012)
%
% Functions
%   - fun_generate_sino_ms
%   - parzenwin_2
%   - fun_interpolation_corr_matrix
%   - fun_gridrec_ms ***(the important one)
%   - fun_ifft_image_corr
%   - fun_calc_error
%
% Note
%   - FFT of a sinogram, size(sino_fft) = [N_layer, N_freq, N_theta]
%   - C is for the interpolation from polar to Cartesian
%   - parzenFilter is applied when filling the Fourier space
%   - Inverse FFT of cartesianGridInterpolatedFFT will be the recon. image
%   - mask is used instead of |k| when filling the Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mex -O COPTIMFLAGS='-O2' LDOPTIMFLAGS='-O2' set_projections_cpu_mex.cpp

clear all; clc; close all; tic

position_range = [0 0 900 800]; FS = 15; 

%% === Input
Np = 256;               % phantom size

N_layer = 13;           % 0 then no multi-slice
use_gridrec = 1;        % 0 (iradon), 1 (both; gridrec slow)
angle_step = 30;        % angular step for tomogram
offset = 0; %0.95       % offset angle array by a little, can ignore
flag_normalize = 1;     % temporary

filter_type = 'Ram-Lak';    % 'Ram-Lak', 'Hann', 'Hamming' , 'parzen' 
C = 7;                      % interpolation (gridrec.py: C=7) *DO NOT CHANGE
widthParam = -1;         % -1 for no filter; For parzenFilter: between 0 and 2, high number means stronger filtering in high freq. domain (gridrec.py: widthParam=1.5)

flag_save_png = 0;      % figure number to save to png
fn_suffix = '_a';       % png filename suffix

param.taper = 20;       % Pixels to taper images - Increase until the FSC does not change anymore
param.SNRt = 0.5;       % SNRt = 0.2071 for 1/2 bit,  0.5 for 1 bit threshold 
param.thickring = 12;
param.align_images = 0;
param.pixelsize = 100; 
param.plot = 2;         % figure number
param2 = param;
param2.plot = 0;

%% === Create phantom
which_phantom = 'mouse';
if strcmp(which_phantom,'sp')
    img_phantom = phantom(Np);
    img_raw = img_phantom;
    zp_len = 10;
    [img_temp, img_orig, fimg_orig] = fun_create_good_phantom(img_phantom, zp_len);
elseif strcmp(which_phantom,'mouse')
    mouse = load('obj2_carbon_3layers_set2.mat');
    ob_dims = size(mouse.obj{2});
    binning = ob_dims(1)/Np;
    img_raw = fftshift(angle(mouse.obj{2}));
    img_raw = img_raw(1:binning:end,1:binning:end);
    zp_len = round(Np/5);   % fairly arbitrary, 5
    [img_phantom, ~, ~] = fun_create_good_phantom(img_raw, zp_len);
    img_orig = img_raw;
    
    ob_dims = size(img_phantom);
    apod = round(zp_len*1.5);   % fairly arbitrary, 1.5
    mask2 = fftshift(fract_hanning_pad(ob_dims(1), ob_dims(1)-zp_len*4+apod, ob_dims(1)-zp_len*4));
    img_phantom = img_phantom.*mask2;
    img_phantom = abs(img_phantom);
    img_phantom = img_phantom./max(img_phantom(:));
    img_temp = img_phantom;
    
    % figure(54); subplot(1,2,1); imagesc(mask2); subplot(1,2,2); imagesc(img_temp); colorbar
end

% === angles

angle_array_true = [0.05+offset:angle_step:179.999+offset];
N_theta = length(angle_array_true);

if numel(angle_array_true)==1 % just because iradonfast_v2 doesn't work for 1 angle
   angle_singleslice = [angle_array_true angle_array_true+0.01];
else
   angle_singleslice = angle_array_true;
end

%% === Generate sino
[sino, pad_dim] = fun_create_good_sino(img_temp, angle_singleslice, Np);
pad_dim = [0 0];
sino = fun_padarray_2D(sino, [pad_dim 0 0], [0 0 0 0]);
N_freq = size(sino,1);

%% === Calculate filters needed
if use_gridrec 
    parzenFilter = parzenwin_2(N_freq, widthParam); 
    alpha = N_theta/3.5; 
%     if N_layer==1
%         alpha = size(sino_fft,3)*N_layer/3.6; %0.99998546;
%     else
%         alpha = pi/(Np/N_layer)*0.15; % no idea why 
%     end
    [lookupTable, Nsupport, tblspcg] = fun_interpolation(C, N_freq);
    interpolationCorrectionMatrix = fun_interpolation_corr_matrix(C, lookupTable, alpha);
end

%% === Single slice
% ---- Using iradonfast_v2 ----
img_recon_iradon = iradonfast_v2(sino, angle_singleslice, filter_type); % Reconstruction
img_recon_iradon_fft = fftshift(fft2(ifftshift(img_recon_iradon)));
[Q_0, img_recon_iradon, ~] = fun_calc_error(img_recon_iradon, img_phantom, param); % correlation, rmse, FSC resolution

% ---- Using Gridrec ----
if use_gridrec
    [sino_gridrec_fft, ~] = fun_generate_sino_ms(angle_array_true, 1, img_temp, pad_dim);
    [img_recon_gridrec_fft, ~] = fun_gridrec_ms(sino_gridrec_fft, angle_array_true, parzenFilter, lookupTable, Nsupport, tblspcg);
    [~,  img_recon_gridrec] = fun_ifft_image_corr(img_recon_gridrec_fft, interpolationCorrectionMatrix);
    if flag_normalize
        max(img_recon_gridrec(:))
        img_recon_gridrec = img_recon_gridrec./max(img_recon_gridrec(:));
    end
    [Q_1, img_recon_gridrec, ~] = fun_calc_error(flipud(img_recon_gridrec), img_phantom, param);
end

% ---- Plot
try
    figure(2)
    subplot(2,3,1);
    % imagesc((abs(img_recon_fine))); colorbar; axis equal tight xy
    % st_title = sprintf(' Fine angles\n theta step %.1f deg \n iradon, filter %s\n (Q, rmse, res) = (%.3f, %.3f, %.3f)', delta, filter_type, Q_fine(1), Q_fine(2), Q_fine(3));
    % title(st_title,'interpreter','none');

    figure(2)
    subplot(2,3,2);
    imagesc((abs(img_recon_iradon))); colorbar; axis equal tight xy
    st_title = sprintf('== iradon ==\n theta step %d deg \n filter %s\n (Q, rmse, res) = (%.3f, %.3f, %.3f)', angle_step, filter_type, Q_0(1), Q_0(2), Q_0(3));
    title(st_title,'interpreter','none');
    if use_gridrec
        subplot(2,3,3);
        imagesc((abs(img_recon_gridrec))); colorbar; axis equal tight xy
        st_title = sprintf('== Gridrec-Single ==\n interpolation C %d \n filter parzen %.1f\n (Q, rmse, res) = (%.3f, %.3f, %.3f)', C, widthParam, Q_1(1), Q_1(2), Q_1(3));
        title(st_title,'interpreter','none');
    end
    set(gcf,'position',position_range,'color','w')
catch
end


%% === Multi-slice: Generate sino and recon
if N_layer>0
    % === Gridrec
    if use_gridrec
        
        angle_array = angle_array_true;
        interpolationCorrectionMatrix = interpolationCorrectionMatrix * (N_theta/length(angle_array)) * Np/N_layer/2;
         
        t0 = tic;
        for loop = [1]
            if loop==1
                [sino_fft, img_ms] = fun_generate_sino_ms(angle_array, N_layer, img_temp, pad_dim);
            else
                [img_now, ~, ~] = fun_create_good_phantom(abs(img_recon), zp_len);
                [angle_array, sino_fft, sino_precession] = fun_generate_sino_ms_precession([0], 0, angle_array, N_layer, img_now, pad_dim);
            end
            
            % === Gridrerc
            [cartesianGridInterpolatedFFT, mask] = fun_gridrec_ms(sino_fft, angle_array, parzenFilter, lookupTable, Nsupport, tblspcg);
            [img_recon,  img_recon_corr] = fun_ifft_image_corr(cartesianGridInterpolatedFFT, interpolationCorrectionMatrix);
            [Q_recon, img_recon, ~] = fun_calc_error(flipud(img_recon), img_phantom, param2);
            [Q_recon_corr, img_recon_corr, img_orig] = fun_calc_error(flipud(img_recon_corr), img_phantom, param);
            if flag_normalize
                temp_max = max(img_recon_corr(:))
                img_recon_corr = img_recon_corr./temp_max;
                cartesianGridInterpolatedFFT = cartesianGridInterpolatedFFT./temp_max;
            end
        end
        t_gridrec = toc(t0);
        
        % ---- Plot results
        figure(2)
        subplot(2,3,6);
        imagesc((abs(img_recon_corr))); colorbar; axis equal tight xy; drawnow
        title(sprintf('---------- N_layer %d ----------\n zp_len = %d, pad_dim = [%d %d]\n (Q, rmse, res) = (%.3f, %.3f, %.3f), %dmin', ...
        N_layer, zp_len, pad_dim(1), pad_dim(2), Q_recon_corr(1), Q_recon_corr(2), Q_recon_corr(3), round(t_gridrec/60)),'interpreter','none');
    end

    %% ---- Plot 
    figure(3);
    img_plot = fun_padarray_2D(img_orig, [pad_dim pad_dim], [0 0 0 0]);
    fimg_temp = fftshift(fft2(ifftshift(img_orig)));
    [img_plot, ~] = fun_crop_images(img_phantom, img_raw);
    subplot(1,3,1); imagesc(img_plot); colorbar; axis equal tight xy; title(sprintf('----- True Image -----\n(before zeropadding)')); 
    caxis([0 1]);
    
    if use_gridrec
        [img_plot, ~] = fun_crop_images(img_recon_gridrec, img_raw);
        subplot(1,3,2); imagesc(img_plot); colorbar; axis equal tight xy; title(sprintf('----- Gridrec-Single -----\n(theta %.1f deg)', angle_step),'interpreter','none');
        caxis([0 1]);
    end
        
    
    [img_plot, ~] = fun_crop_images(img_recon_corr, img_raw);
    subplot(1,3,3); imagesc(img_plot); colorbar; axis equal tight xy; title(sprintf('----- Gridrec-MS -----\n(N_layer %d) ', N_layer),'interpreter','none');
    caxis([0 1]); colormap gray
    set(gcf,'position',position_range,'color','w')

    
    figure(4);
    subplot(1,3,1); imagesc(log10(abs(fimg_temp))); colorbar; axis equal tight xy; 
    title(sprintf('----- log10(abs(FFT(True Image)) -----\n')); 
    caxis([0 4]);
    plot_range = caxis;  
    
    fimg_temp = fftshift(fft2(ifftshift(img_recon_gridrec)));
    subplot(1,3,2); imagesc(log10(abs(fimg_temp))); colorbar; axis equal tight xy; 
    title(sprintf('----- Gridrec-Single -----\n(theta %.1f deg)', angle_step),'interpreter','none');
    caxis(plot_range);

    [img_plot, ~] = fun_crop_images(img_recon_corr, img_raw);
    fimg_temp = fftshift(fft2(ifftshift(img_recon_corr)));
    subplot(1,3,3); imagesc(log10(abs(fimg_temp)));  colorbar; axis equal tight xy; 
    title(sprintf('----- log10(abs(FFT)), Gridrec-MS -----\n(N_layer %d) ', N_layer),'interpreter','none');
    caxis(plot_range);
    
    colormap franzmap
    set(gcf,'position',position_range,'color','w')
    
end

%% === Plot 
f1 = figure(1); clf; M=4; N=6; 

% =======
fig = 1;
subplot(M,N,fig);
imagesc((abs(img_recon_iradon))); colorbar; axis equal tight xy
st_title = sprintf('== iradonfast_v2 ==\n single slice\n (Q, rmse, res) = (%.3f, %.3f, %.3f)', Q_0(1), Q_0(2), Q_0(3));
title(st_title,'interpreter','none');
%
subplot(M,N,fig+N);
imagesc((angle(img_recon_iradon))); colorbar; axis equal tight xy
title('phase','interpreter','none');
%
subplot(M,N,fig+N*2);
imagesc(log10(abs(img_recon_iradon_fft))); colorbar; axis equal tight xy
title('FFT(recon image).abs','interpreter','none');
%
subplot(M,N,fig+N*3);
imagesc((angle(img_recon_iradon_fft))); colorbar; axis equal tight xy
title('phase','interpreter','none');

% =======
if use_gridrec
    fig = 2;
    subplot(M,N,fig);
    imagesc((abs(img_recon_gridrec))); colorbar; axis equal tight xy
    st_title = sprintf('== Gridrec ==\n single slice\n (Q, rmse, res) = (%.3f, %.3f, %.3f)', Q_1(1), Q_1(2), Q_1(3));
    title(st_title,'interpreter','none');
    %
    subplot(M,N,fig+N);
    imagesc((angle(img_recon_gridrec))); colorbar; axis equal tight xy
    title('phase','interpreter','none');
    %
    subplot(M,N,fig+N*2);
    imagesc(log10(abs(img_recon_gridrec_fft))); colorbar; axis equal tight xy
    title('FFT(recon image).abs','interpreter','none');
    %
    subplot(M,N,fig+N*3);
    imagesc((angle(img_recon_gridrec_fft))); colorbar; axis equal tight xy
    title('phase','interpreter','none');
end

% =======
fig = 3;
subplot(M,N,fig);
imagesc(img_temp); colorbar; axis equal tight xy
st_title = sprintf('True Image\nNp=%d\nN_layer=%d', Np, N_layer);
title(st_title,'interpreter','none','fontsize',15);

% =======
if use_gridrec
    fig = 4;
    subplot(M,N,fig);
        plot(lookupTable,'b.'); grid on; axis tight
        st_title = sprintf('(PSWF) Convolvent\nlength = %d',length(lookupTable));
        title(st_title,'interpreter','none'); 
    %
    subplot(M,N,fig+N);
    for q = 1:2:N_freq
        for angle_temp = [0:30:180]
            k = q - N_freq/2;
            theta = angle_temp/180*pi;

            kx_p = k*cos(theta);
            ky_p = k*sin(theta);
            hold on; plot(kx_p, ky_p,'r.'); grid on;
        end
    end
    grid on; xlabel('kx'); ylabel('ky'); box on; axis tight equal
    title(sprintf('Map sino_fft \npolarToCartesian \n e.g. [0, 45, 90] deg'),'interpreter','none'); 
    %
    subplot(M,N,fig+N*2); 
    for q = 1:round(N_freq/2)
        for angle_temp = angle_array_true(round(N_theta/2))
            k = q - N_freq/2;
            theta = angle_temp/180*pi;

            % polarToCartesian
            kx_p = k*cos(theta);
            ky_p = k*sin(theta);

            kx_p_NN = round(kx_p);
            ky_p_NN = round(ky_p);

            local_kx_range_min = kx_p_NN - Nsupport;
            local_ky_range_min = ky_p_NN - Nsupport;
            local_kx_range_max = kx_p_NN + Nsupport;% + 1;
            local_ky_range_max = ky_p_NN + Nsupport;% + 1;

            kx = [local_kx_range_min:local_kx_range_max];
            ky = [local_ky_range_min:local_ky_range_max];

            delta_kx = kx_p - kx;
            delta_ky = ky_p - ky;
            convolveArg_x = tblspcg*abs(delta_kx);
            convolveArg_y = tblspcg*abs(delta_ky);

            weight_kx = lookupTable(round(convolveArg_x)+1); 
            weight_ky = lookupTable(round(convolveArg_y)+1);
            weight = weight_ky'*weight_kx; 

        end
    end
    im = imagesc(kx,ky,weight); colorbar; axis equal tight xy 
    %im.AlphaData = 0.5;
    hold on; plot(kx_p, ky_p,'rx','markersize',6,'linewidth',2);
    grid on; xlabel('kx'); ylabel('ky'); axis equal tight xy 
    [KX, KY] = meshgrid(kx,ky);
    hold on; plot(KX,KY,'wo','markersize',2,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor','w');
    title(sprintf('weight, C = %d\ncontribution = \nsino_fft(q,theta)*weight*|k|',C),'interpreter','none');
    %
    subplot(M,N,fig+N*3);
    q = 1:N_freq;
    k = q - N_freq/2;
    plot(parzenFilter,'.','color',[0.1 0.6 0.2]); grid on; axis tight
    %hold on; plot(parzenFilter(:).*abs(k(:)),'color',[0.1 0.8 0.2]); axis tight
    st_title = sprintf('parzenFilter\nwidthParam = %.1f, length = %d\ncontribution*parzenFilter( round(|q|) )',widthParam,length(parzenFilter));
    title(st_title,'interpreter','none');

    % =======
    fig = 5;
    subplot(M,N,fig);
    imagesc((abs(img_recon))); colorbar; axis equal tight xy
    %title(sprintf('img_recon \n (Q, rmse, res) = (%.3f, %.3f, %.3f)', Q_recon(1), Q_recon(2), Q_recon(3)),'interpreter','none');
    subplot(M,N,fig+N);
    imagesc((angle(img_recon))); colorbar; axis equal tight xy
    title('phase','interpreter','none');
    subplot(M,N,fig+N*2);
    imagesc(log10(abs(cartesianGridInterpolatedFFT))); colorbar; axis equal tight xy
    st_title = sprintf('FFT\nsize=(%d, %d)',size(cartesianGridInterpolatedFFT,1), size(cartesianGridInterpolatedFFT,2));
    title(st_title,'interpreter','none');
    subplot(M,N,fig+N*3);
    imagesc((angle(cartesianGridInterpolatedFFT))); colorbar; axis equal tight xy
    title('phase','interpreter','none');
end

% =======
fig = 6;
subplot(M,N,fig);
imagesc((abs(img_recon_corr))); colorbar; axis equal tight xy; %caxis([0 500])
if use_gridrec==0
    title(sprintf('== irandon ==\nimg_recon_corr \n (Q, rmse, res) = (%.3f, %.3f, %.3f)', Q_recon_corr, Q_recon_corr(2), Q_recon_corr(3)),'interpreter','none');
else
    title(sprintf('== GridrecMS ==\nimg_recon_corr \n (Q, rmse, res) = (%.3f, %.3f, %.3f)', Q_recon_corr, Q_recon_corr(2), Q_recon_corr(3)),'interpreter','none');
end
if use_gridrec
    subplot(M,N,fig+N);
    surfc(((interpolationCorrectionMatrix)),'Edgecolor','none'); colorbar; %axis equal tight xy
    title(sprintf('interpolationCorrectionMatrix \n (pi/2/C/alpha) = %.1f',pi/2/C/alpha),'interpreter','none');

    subplot(M,N,fig+N*2);
    imagesc(mask); colorbar; axis equal tight xy
    title('mask','interpreter','none');
end

% =======
set(gcf,'position',[0 50 1800 1000],'color','w')

toc



%% === Save to png
figure(2); colormap gray
if flag_save_png>0
    tic
    [status, loc]= system('pwd');          
    if flag_save_png==1
        filename = sprintf('%s/figs/tomo_syn_Np%d_Nlayer%d_%.0fdeg_delta%0.2f_%s_fig1.png',loc(1:end-1), Np, N_layer, angle_step, delta, fn_suffix);
        print('-f1','-dpng','-r300',filename);
    elseif flag_save_png==2
        filename = sprintf('%s/figs/tomo_syn_Np%d_Nlayer%d_%.0fdeg_delta%0.2f_%s_fig2_gray.png',loc(1:end-1), Np, N_layer, angle_step, delta, fn_suffix);
        print('-f2','-dpng','-r300',filename);
    end
    fprintf('saving %s\n',filename);   
    toc
end


