%% Step 2a: Tune Bg_Sigma before running the full correlation-image pipeline
% Run this once at the start of your analysis (or once per animal, if
% zoom/resolution differs between recordings) to pick a good Bg_Sigma for
% compute_correlation_image.m / Step2_CorrelationImage_and_ROIs.m.

Repo_Path = 'C:\Users\afermigier\Documents\GitHub\1Photon_Analysis';
addpath(genpath(Repo_Path));

mouse_id = '840';   % without the G/R suffix
task     = 'EPM';
base_dir = 'F:\Inscopix_Projects\202508_DualColorMiniscope';

green_mouse_color = [mouse_id 'G'];
red_mouse_color   = [mouse_id 'R'];
green_channel_label = [green_mouse_color '_' task];
red_channel_label   = [red_mouse_color '_' task];

green_out_folder = fullfile(base_dir, green_mouse_color, green_channel_label, 'processed_data');
red_out_folder   = fullfile(base_dir, red_mouse_color,   red_channel_label,   'processed_data');

% Green's mean-padded output is MC_meanpad.mat (kept separate from the
% pipeline's own zero-padded MC.mat); red's is saved as MC.mat itself.
green_mc_mat_path = fullfile(green_out_folder, 'MC_meanpad.mat');
red_mc_mat_path   = fullfile(red_out_folder,   'MC.mat');

Frame_Start = 200;   % a bit into the recording, same reasoning as the QC snapshots
N_Frames = 600;      % enough for a stable estimate; this is fast either way

%% 1) Measure a typical neuron radius on the green channel (best-resolved cells)
MC_f_G = matfile(green_mc_mat_path);
green_preview = mean(double(MC_f_G.MC(:, :, Frame_Start:Frame_Start+199)), 3);

radius_px = measure_neuron_radius(green_preview, prctile(green_preview(:), [1 99]));
fprintf('\n==> Measured neuron radius: %.1f px\n', radius_px);
fprintf('==> Starting point for Bg_Sigma: roughly %.1f - %.1f px (1-1.5x the radius)\n\n', ...
        radius_px, 1.5*radius_px);

%% 2) Sweep Bg_Sigma around that estimate and compare visually - green channel
Bg_Sigma_List = round(radius_px * [0.5 0.75 1 1.5 2]);
fprintf('Previewing Bg_Sigma values on GREEN: %s\n', mat2str(Bg_Sigma_List));
preview_correlation_image(green_mc_mat_path, Frame_Start, N_Frames, Bg_Sigma_List);
sgtitle_txt = 'Green channel - Bg\_Sigma sweep';
annotation(gcf, 'textbox', [0 0.96 1 0.03], 'String', sgtitle_txt, 'Interpreter', 'tex', ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontWeight', 'bold');

%% 3) Same sweep on the red channel - the noisier channel is the real test
fprintf('Previewing Bg_Sigma values on RED: %s\n', mat2str(Bg_Sigma_List));
preview_correlation_image(red_mc_mat_path, Frame_Start, N_Frames, Bg_Sigma_List);
annotation(gcf, 'textbox', [0 0.96 1 0.03], 'String', 'Red channel - Bg\_Sigma sweep', ...
    'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontWeight', 'bold');

% Pick the smallest Bg_Sigma where neurons look like clean, well-separated,
% roughly round blobs on BOTH channels (they should use the same value -
% it reflects the optics/resolution, not the channel). Too small: cells
% look dim/fragmented (their own signal leaking into the background
% estimate). Too large: blobs merge, or broad illumination structure
% starts to reappear. Set that chosen value as Bg_Sigma in
% Step2_CorrelationImage_and_ROIs.m.