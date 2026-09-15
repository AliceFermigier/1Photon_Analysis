%% Single-session test: apply green MC shifts to the red channel
% Edit the paths below to match your 839G_EPM / 839R_EPM session, then run.

Repo_Path = 'C:\Users\afermigier\Documents\GitHub\1Photon_Analysis';
addpath(genpath(Repo_Path));

Output_Path = 'C:\Users\afermigier\Documents\GitHub\1Photon_Analysis\Functions';
addpath(Output_Path);

green_tif  = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data\840G\840G_RewardAirpuff\840G_RewardAirpuff.tiff';
red_tif    = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data\840R\840R_RewardAirpuff\840R_RewardAirpuff.tiff';
shifts_mat = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data\840G\840G_RewardAirpuff\processed_data\MC_Shifts.mat';
out_folder = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data\840R\840R_RewardAirpuff';

Chunk_size = 3000;
Interleave_Mode = 'auto';

apply_shifts_to_red(green_tif, shifts_mat, red_tif, out_folder, Chunk_size, Interleave_Mode);

%% Quick visual sanity check
MC_f = matfile(fullfile(out_folder, 'MC.mat'));
raw_R = loadtiff(red_tif, 1, 200);   % first 200 raw red frames

figure;
subplot(1,2,1); imagesc(mean(raw_R, 3)); axis image off; colormap gray;
title('Raw red - mean of first 200 frames');

subplot(1,2,2); imagesc(mean(MC_f.MC(:,:,1:200), 3)); axis image off; colormap gray;
title('Motion-corrected red - mean of first 200 frames');

% If motion correction worked, the corrected mean image should look
% noticeably sharper / less blurred than the raw one, since raw motion
% blurs the mean image over time.

% Load one representative session
MC_f = matfile(fullfile(out_folder, 'MC.mat'));
raw_R = loadtiff(red_tif, 1, 200);
corr_R = MC_f.MC(:,:,1:200);

% --- 1) Compare a SINGLE frame, not an average, with matched color scale ---
f = 100;
figure;
subplot(1,2,1); imagesc(raw_R(:,:,f));  axis image off; colormap gray; title('raw, frame 100');
subplot(1,2,2); imagesc(corr_R(:,:,f)); axis image off; colormap gray; title('corrected, frame 100');
clim_common = [min(raw_R(:)) max(raw_R(:))];
subplot(1,2,1); caxis(clim_common);
subplot(1,2,2); caxis(clim_common);   % same scale on both now

% --- 2) Basic stats ---
fprintf('Raw:  min=%d max=%d mean=%.1f std=%.1f\n', min(raw_R(:)), max(raw_R(:)), mean(raw_R(:)), std(double(raw_R(:))));
fprintf('Corr: min=%d max=%d mean=%.1f std=%.1f\n', min(corr_R(:)), max(corr_R(:)), mean(corr_R(:)), std(double(corr_R(:))));

% --- 3) Standard-deviation projection - this is the real test.
% If registration is working, cells should appear as sharper, higher-contrast
% blobs in the corrected std-projection than the raw one (motion smears a
% moving cell's std signal across its trajectory; a stationary cell has a
% tight, high std footprint right on the soma).
figure;
subplot(1,2,1); imagesc(std(double(raw_R), 0, 3));  axis image off; title('raw std projection');
subplot(1,2,2); imagesc(std(double(corr_R), 0, 3)); axis image off; title('corrected std projection');