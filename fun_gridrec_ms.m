function [cartesianGridInterpolatedFFT, mask_use] = fun_gridrec_ms(sino_fft, angle_array, parzenFilter, lookupTableOfConvolventInFourierSpace, numberOfSupportNeighbours, tblspcg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2017-01-19
%
% Original: see gridRec.py from TOMCAT and 
% F. Marone and M. Stampanoni, “Regridding reconstruction algorithm for 
% real-time tomographic imaging,�? J.Synchrotron Radiat.19, 1029–1037 (2012)
%
% Input
%   - FFT of a sinogram, size(sino_fft) = [N_layer, N_freq, N_theta]
%   - C is for the interpolation from polar to Cartesian
%   - parzenFilter is applied when filling the Fourier space
% Output
%   - Inverse FFT of cartesianGridInterpolatedFFT will be the recon. image
%   - mask is used instead of |k| when filling the Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_use_interp1 = 0;       % can use 'interp1' instead of 'round', but it's slow
flag_plot_mask_anime = 0;   % flag, as well as figure number

N_theta = length(angle_array);
 
N_layer = size(sino_fft,1);
N_freq = size(sino_fft,2);
bw_freq = ceil((N_freq+1)/2); % 7-->4; 8-->5

cartesianGridInterpolatedFFT = complex(zeros(N_freq,N_freq, 'single'));
mask = complex(zeros(N_freq,N_freq, 'single')); %zeros(N_freq,N_freq);

%[lookupTableOfConvolventInFourierSpace, numberOfSupportNeighbours, tblspcg] = fun_interpolation(C, N_freq);

%% Gridrec
%idx=1
fprintf('----- Gridrec (N_layer = %d) -----\n', N_layer);
for thetaIndex = 1:N_theta
    theta = angle_array(thetaIndex)/180*pi;

    if mod(thetaIndex,30)==0; fprintf('%.0f %%\n',thetaIndex/N_theta*100); end
    %---- Testing different mask
    cartesianGridInterpolatedFFT_old_old = cartesianGridInterpolatedFFT;
    if mod(round(angle_array(thetaIndex)),90)==0
        coeff = 1.0;
        fprintf('theta %.2f deg, coeff %.2f\n',angle_array(thetaIndex), coeff);
    else
        coeff = 1.0;
    end
    %----
    
    for q = 1:N_freq % index for frequency        
        for n =  [1:N_layer]
            k = q - bw_freq; 
            
            % polarToCartesian
            kx_p = k*cos(theta);
            ky_p = k*sin(theta);
            
            distance = n - (N_layer+1)/2;
            kx_p = kx_p - distance*cos(pi/2-theta);
            ky_p = ky_p + distance*sin(pi/2-theta);
            k_corr = sqrt(k^2 + distance^2); %?

            kx_p_NN = round(kx_p);
            ky_p_NN = round(ky_p);

            % To check the if conditions is not necessary anymore, since ... - 2*numberOfSupportNeighbours in ref (1)
            if kx_p_NN^2 + ky_p_NN^2 > (N_freq/2 - 2*numberOfSupportNeighbours)^2
                continue
            else
                %fprintf('stopped\n');
            end
                               
            local_kx_range_min = kx_p_NN - numberOfSupportNeighbours;
            local_ky_range_min = ky_p_NN - numberOfSupportNeighbours;
            local_kx_range_max = kx_p_NN + numberOfSupportNeighbours;% + 1;
            local_ky_range_max = ky_p_NN + numberOfSupportNeighbours;% + 1;

            kx = [local_kx_range_min:local_kx_range_max];
            ky = [local_ky_range_min:local_ky_range_max];
            % calculate the distance between a potentially contributing point and the cartesian sampling point (kx, ky)
            delta_kx = kx_p - kx;
            delta_ky = ky_p - ky;
            convolveArg_x = tblspcg*abs(delta_kx);
            convolveArg_y = tblspcg*abs(delta_ky);
            % calculate weight arrays and matrix
            if flag_use_interp1
                weight_kx = interp1(1:length(lookupTableOfConvolventInFourierSpace), lookupTableOfConvolventInFourierSpace, (convolveArg_x)+1);
            	weight_ky = interp1(1:length(lookupTableOfConvolventInFourierSpace), lookupTableOfConvolventInFourierSpace, (convolveArg_y)+1);
            else
                weight_kx = lookupTableOfConvolventInFourierSpace(round(convolveArg_x)+1);
                weight_ky = lookupTableOfConvolventInFourierSpace(round(convolveArg_y)+1);
            end
            weight = weight_ky'*weight_kx; % convolution 2D-mask (tensor product of 1D convolution arrays, psi_fft(x,y) = psi_fft(x) psi_fft(y) )
            
            % multiplication with ramp filter numpy.abs(k) is an idealization and does not take discrete 
            % nature of FFT into account (See regridding paper by F. Marone)
            
            if 1 
                if N_layer==1 
                    contribution = sino_fft(n, q, thetaIndex)*weight*abs(k);            
                else
                    if distance ==0
                        dd = 1;
                    else
                        dd = abs(distance);
                    end
                    contribution = sino_fft(n, q, thetaIndex)*weight*abs(1);     
                end
            else
                contribution = zeros(size(weight)); 
            end
            
            if flag_use_interp1
                temp_parzen = interp1(1:length(parzenFilter),  parzenFilter, abs(k_corr+ceil((N_freq+1)/2)));
            else
                temp_parzen = parzenFilter( round( abs(k_corr+ceil((N_freq+1)/2)) ) );
            end
            %filteredContribution(:,:,idx) = single(contribution*temp_parzen);
            filteredContribution = contribution*temp_parzen;
            
            rowIndex = ceil(ky + N_freq/2)+0;   % .astype(numpy.int) is floor in python
            columnIndex = ceil(kx + N_freq/2)+0;  
            %positions(idx,:) = [rowIndex(1)-1, columnIndex(1)-1];
            cartesianGridInterpolatedFFT_old = cartesianGridInterpolatedFFT;
            cartesianGridInterpolatedFFT(rowIndex, columnIndex) ...
                = cartesianGridInterpolatedFFT(rowIndex, columnIndex) + filteredContribution;
            %mask_weight(rowIndex, columnIndex) = mask_weight(rowIndex, columnIndex) + weight;
            
            if N_layer>1 
                mask = mask + coeff*(abs(cartesianGridInterpolatedFFT-cartesianGridInterpolatedFFT_old) > 1e-9); % increase val (now 0) will enhance the high frequencies      
                if 0 %flag_plot_mask_anime
                    figure(flag_plot_mask_anime); imagesc(mask); drawnow; 
                    %keyboard
                end
            end
            
            %idx = idx+1

        end
         
    end
    
    if flag_plot_mask_anime && N_layer>1 
        figure(flag_plot_mask_anime); clf
        imagesc(mask);   
        axis tight equal xy; colorbar; caxis([0 138])
        set(gca,'fontsize',20,'fontweight','bold')
        
        figure(flag_plot_mask_anime+1); 
        imagesc(log10(abs(cartesianGridInterpolatedFFT)));
        axis tight equal xy; colorbar; caxis([-2 4])
        set(gca,'fontsize',20,'fontweight','bold')
        drawnow; 
        %keyboard; 
    end
end
   
%%% If using mex 
% positions = positions(1:size(filteredContribution,3),:);
% set_projections_cpu_mex(cartesianGridInterpolatedFFT, filteredContribution, int32(positions));
% temp_mask = abs(filteredContribution);
% temp_mask((temp_mask>0))=1;
% temp_mask((temp_mask<=0))=0;
% set_projections_cpu_mex(mask, temp_mask, int32(positions));  


if N_layer > 1 
    
    mask = real(mask); 
    figure(100+N_layer); imagesc(mask); colorbar
    
    mask(mask==0) = 1;
    mask = mask./max(mask(:));
    

    mask_use = mask; %%
    cartesianGridInterpolatedFFT = cartesianGridInterpolatedFFT./mask_use;

    figure(99); 
    imagesc(log10(abs(cartesianGridInterpolatedFFT)));
    axis tight equal xy; colorbar; caxis([-2 4])
    set(gca,'fontsize',20,'fontweight','bold')
    drawnow; 

else
    mask_use = ones(size(mask));
end

