%% Step 2: Correlation-image background, ROI drawing, and fluorescence extraction
% Run this once per session per channel (edit Channel/paths below).
%
% Supersedes Step2_MaxProjection_and_ROIs.m: the background for ROI
% drawing is now a local pixel-correlation image (same underlying concept
% CNMFE uses for its own Outlines_Neurons.png), not a max-intensity
% projection. See compute_correlation_image.m for why this makes neurons
% stand out more clearly and why it doesn't need a separate global
% illumination correction (that correction, and compute_max_projection.m,
% are still available standalone if you ever want a plain intensity view).

Repo_Path = 'C:\Users\afermigier\Documents\GitHub\1Photon_Analysis';
addpath(genpath(Repo_Path));

mouse_color = '840R';
task        = 'EPM';

% --- Pick the channel/session to work on ---
base_dir      = 'F:\Inscopix_Projects\202508_DualColorMiniscope';
channel_label = [mouse_color '_' task];
out_folder    = fullfile(base_dir, mouse_color, channel_label, 'processed_data');

% Red's own mean-padded output from step 1 is saved as MC.mat; green's is
% saved separately as MC_meanpad.mat so it never overwrites the pipeline's
% own zero-padded MC.mat (see apply_shifts_to_red.m / apply_shifts_to_green.m).
if contains(mouse_color, 'R')
    mc_mat_path         = fullfile(out_folder, 'MC.mat');
    shifts_applied_path = fullfile(out_folder, 'MC_Shifts_applied.mat');
else
    mc_mat_path         = fullfile(out_folder, 'MC_meanpad.mat');
    shifts_applied_path = fullfile(out_folder, 'MC_meanpad_Shifts_applied.mat');
end

Chunk_size = 3000;
Bg_Sigma = 15;   % pixels - set a bit larger than a neuron's radius (see
                 % compute_correlation_image.m for what this controls)

%% Compute the local correlation image
Cn = compute_correlation_image(mc_mat_path, Chunk_size, Bg_Sigma);

save(fullfile(out_folder, [channel_label '_CorrImg.mat']), 'Cn', 'Bg_Sigma');

% Quick 8-bit PNG for viewing outside MATLAB (Cn is a correlation map,
% roughly in [-1, 1], so it's rescaled to [0, 255] rather than saved as a
% raw-intensity-style 16-bit tif).
Cn_disp = Cn - min(Cn(:));
Cn_disp = uint8(255 * Cn_disp / max(Cn_disp(:)));
imwrite(Cn_disp, fullfile(out_folder, [channel_label '_CorrImg.png']));

%% Visual check
figure('Name', 'Local correlation image');
imagesc(Cn, prctile(Cn(:), [1 99.5])); axis image off; colormap gray; colorbar;
title(sprintf('%s - local correlation image (Bg\\_Sigma = %d)', channel_label, Bg_Sigma), ...
      'Interpreter', 'tex');

%% Draw ROIs on the correlation image
clim = prctile(Cn(:), [1 99.5]);
ROIs = draw_polygon_rois(Cn, clim);

save(fullfile(out_folder, [channel_label '_ROIs.mat']), 'ROIs');

%% Extract fluorescence traces (from the RAW, uncorrected MC data) and save to CSV
out_csv = fullfile(out_folder, [channel_label '_Traces.csv']);
Traces = extract_roi_traces(mc_mat_path, ROIs, Chunk_size, out_csv);

%% Quick check: plot a couple of extracted traces
figure('Name', 'Example extracted traces');
n_show = min(4, numel(ROIs));
for k = 1:n_show
    subplot(n_show, 1, k);
    plot(Traces(:, k));
    ylabel(sprintf('Neuron %d', ROIs(k).ID));
    if k == 1, title('Example fluorescence traces (raw, mean pixel value per ROI)'); end
end
xlabel('Frame');