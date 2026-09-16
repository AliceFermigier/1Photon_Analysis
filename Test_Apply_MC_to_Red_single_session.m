%% Single-session test: apply green MC shifts to the red channel

Repo_Path = 'C:\Users\afermigier\Documents\GitHub\1Photon_Analysis';
addpath(genpath(Repo_Path));

Output_Path = 'C:\Users\afermigier\Documents\GitHub\1Photon_Analysis\Functions';
addpath(Output_Path);

green_tif  = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data\840G\840G_RewardAirpuff\840G_RewardAirpuff.tiff';
red_tif    = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data\840R\840R_RewardAirpuff\840R_RewardAirpuff.tiff';
shifts_mat = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data\840G\840G_RewardAirpuff\processed_data\MC_Shifts.mat';
out_folder = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data\840R\840R_RewardAirpuff';
green_out_folder = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data\840G\840G_RewardAirpuff\processed_data';

Chunk_size = 3000;
Interleave_Mode = 'auto';

%%
apply_shifts_to_red(green_tif, shifts_mat, red_tif, out_folder, Chunk_size, Interleave_Mode);

%% Quick visual sanity check
MC_f = matfile(fullfile(out_folder, 'MC.mat'));
corr_R = MC_f.MC(:,:,1:200);
raw_R = loadtiff(red_tif, 1, 200);   % first 200 raw red frames

MC_f_G = matfile(fullfile(green_out_folder, 'MC.mat'));
corr_G = MC_f_G.MC(:,:,1:200);
raw_G = loadtiff(green_tif, 1, 200);

figure;
subplot(1,2,1); imagesc(mean(raw_R, 3)); axis image off; colormap gray;
title('Raw red - mean of first 200 frames');

subplot(1,2,2); imagesc(mean(corr_R, 3)); axis image off; colormap gray;
title('Motion-corrected red - mean of first 200 frames');

figure;
subplot(1,2,1); imagesc(mean(raw_G, 3)); axis image off; colormap gray;
title('Raw green - mean of first 200 frames');

subplot(1,2,2); imagesc(mean(corr_G, 3)); axis image off; colormap gray;
title('Motion-corrected green - mean of first 200 frames');

% If motion correction worked, the corrected mean image should look
% noticeably sharper / less blurred than the raw one, since raw motion
% blurs the mean image over time.

% --- Basic stats ---
fprintf('Raw:  min=%d max=%d mean=%.1f std=%.1f\n', min(raw_R(:)), max(raw_R(:)), mean(raw_R(:)), std(double(raw_R(:))));
fprintf('Corr: min=%d max=%d mean=%.1f std=%.1f\n', min(corr_R(:)), max(corr_R(:)), mean(corr_R(:)), std(double(corr_R(:))));

% --- Standard-deviation projection - this is the real test.
% If registration is working, cells should appear as sharper, higher-contrast
% blobs in the corrected std-projection than the raw one (motion smears a
% moving cell's std signal across its trajectory; a stationary cell has a
% tight, high std footprint right on the soma).
figure;
subplot(1,2,1); imagesc(std(double(raw_R), 0, 3));  axis image off; title('Raw red std projection');
subplot(1,2,2); imagesc(std(double(corr_R), 0, 3)); axis image off; title('Corrected red std projection');

figure;
subplot(1,2,1); imagesc(std(double(raw_G), 0, 3));  axis image off; title('Raw green std projection');
subplot(1,2,2); imagesc(std(double(corr_G), 0, 3)); axis image off; title('Corrected green std projection');

%% Diagnose the "washed out" corrected std projection
% Run this right after loading raw_R and corr_R (or MC_f.MC) as before.

std_raw  = std(double(raw_R), 0, 3);
std_corr = std(double(corr_R), 0, 3);



%% 1) Look at the actual numeric range of both std maps
fprintf('std_raw:  min=%.2f max=%.2f median=%.2f\n', min(std_raw(:)),  max(std_raw(:)),  median(std_raw(:)));
fprintf('std_corr: min=%.2f max=%.2f median=%.2f\n', min(std_corr(:)), max(std_corr(:)), median(std_corr(:)));

% If std_corr's max is dramatically higher than its median/typical value,
% a few extreme pixels (almost certainly the flickering zero-padded
% border) are stretching the color scale and crushing everything else
% toward the bottom of the colormap.

%% 2) Where are the extreme values? (should be a thin ring around the edges)
[~, idx] = max(std_corr(:));
[r, c] = ind2sub(size(std_corr), idx);
fprintf('Max std_corr pixel is at row=%d, col=%d (image is %d x %d)\n', r, c, size(std_corr,1), size(std_corr,2));

%% 3) Re-plot with a border crop and a percentile-based color scale
margin = 15;   % pixels to crop from each edge - adjust based on your max shift magnitude
H = size(std_corr,1); W = size(std_corr,2);
crop_rows = (margin+1):(H-margin);
crop_cols = (margin+1):(W-margin);

std_raw_crop  = std_raw(crop_rows, crop_cols);
std_corr_crop = std_corr(crop_rows, crop_cols);

clim_lo = min(prctile(std_raw_crop(:),1), prctile(std_corr_crop(:),1));
clim_hi = max(prctile(std_raw_crop(:),99), prctile(std_corr_crop(:),99));

figure;
subplot(1,2,1); imagesc(std_raw_crop, [clim_lo clim_hi]);  axis image off; colorbar;
title('raw std, cropped, shared scale');
subplot(1,2,2); imagesc(std_corr_crop, [clim_lo clim_hi]); axis image off; colorbar;
title('corrected std, cropped, shared scale');

% Now: do you see cell-shaped high-std blobs in the corrected panel,
% ideally tighter/more punctate than in the raw panel? That would confirm
% registration is genuinely working and the earlier "all blue" look was
% just the border artifact stretching the color scale.