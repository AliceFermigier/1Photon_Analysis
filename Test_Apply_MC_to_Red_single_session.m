%% Single-session test: apply green MC shifts to the red channel
% Edit the paths below to match your 839G_EPM / 839R_EPM session, then run.

Repo_Path = 'C:\Users\afermigier\Documents\GitHub\1Photon_Analysis';
addpath(genpath(Repo_Path));

Output_Path = 'C:\Users\afermigier\Documents\GitHub\1Photon_Analysis\Functions';
addpath(Output_Path);

green_tif  = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data\840G\840G_RewardAirpuff\840G_RewardAirpuff.tiff';
red_tif    = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data\840R\840R_RewardAirpuff\840R_RewardAirpuff.tiff';
shifts_mat = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data\840G\840G_RewardAirpuff\processed_data\MC_Shifts.mat';
out_folder = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data\840R\840R_RewardAirpuff\processed_data';
green_out_folder = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data\840G\840G_RewardAirpuff\processed_data';

Chunk_size = 3000;
Interleave_Mode = 'auto';
Fill_Value = 'frame_mean';   % 'frame_mean' (recommended) | 'zero' | numeric constant

% --- Red channel: shifts borrowed from green, interleave-matched ---
apply_shifts_to_red(green_tif, shifts_mat, red_tif, out_folder, Chunk_size, Interleave_Mode, Fill_Value);

% --- Green channel: regenerated with the SAME Fill_Value, for a fair comparison.
% Saved as MC_meanpad.mat - does NOT touch the pipeline's original MC.mat.
apply_shifts_to_green(green_tif, shifts_mat, green_out_folder, Chunk_size, Fill_Value);

%% Quick visual sanity check - red
MC_f = matfile(fullfile(out_folder, 'MC.mat'));
raw_R  = loadtiff(red_tif, 1, 200);   % first 200 raw red frames
corr_R = MC_f.MC(:,:,1:200);          % first 200 corrected red frames

figure('Name', 'Red channel');
subplot(1,2,1); imagesc(mean(raw_R, 3)); axis image off; colormap gray;
title('Raw red - mean of first 200 frames');
subplot(1,2,2); imagesc(mean(corr_R, 3)); axis image off; colormap gray;
title('Motion-corrected red - mean of first 200 frames');

fprintf('Raw R:  min=%d max=%d mean=%.1f std=%.1f\n', min(raw_R(:)), max(raw_R(:)), mean(raw_R(:)), std(double(raw_R(:))));
fprintf('Corr R: min=%d max=%d mean=%.1f std=%.1f\n', min(corr_R(:)), max(corr_R(:)), mean(corr_R(:)), std(double(corr_R(:))));

%% Quick visual sanity check - green (same Fill_Value as red, for comparison)
MC_f_G = matfile(fullfile(green_out_folder, 'MC_meanpad.mat'));
raw_G  = loadtiff(green_tif, 1, 200);   % first 200 raw green frames
corr_G = MC_f_G.MC(:,:,1:200);          % first 200 corrected green frames

figure('Name', 'Green channel');
subplot(1,2,1); imagesc(mean(raw_G, 3)); axis image off; colormap gray;
title('Raw green - mean of first 200 frames');
subplot(1,2,2); imagesc(mean(corr_G, 3)); axis image off; colormap gray;
title('Motion-corrected green (mean-pad) - mean of first 200 frames');

fprintf('Raw G:  min=%d max=%d mean=%.1f std=%.1f\n', min(raw_G(:)), max(raw_G(:)), mean(raw_G(:)), std(double(raw_G(:))));
fprintf('Corr G: min=%d max=%d mean=%.1f std=%.1f\n', min(corr_G(:)), max(corr_G(:)), mean(corr_G(:)), std(double(corr_G(:))));

% For green specifically, the corrected image should look visibly sharper
% than raw - green has clear, well-visible pyramidal neurons, so this is
% the best channel to visually confirm the whole registration approach
% (interpolation, fill value, chunking) is behaving as expected before
% trusting the same machinery on the noisier red channel.