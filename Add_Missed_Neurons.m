% Run this AFTER Functions/CNMFE.m has produced Result_CNMFE.mat (or
% Result_CNMFE_postprocessed.mat) for a session - works even if CNMFE
% found zero neurons for that session (see post_process_CNMFE.m: that
% Dummyvar/skip branch only triggers in POST-processing; Result_CNMFE.mat
% itself is saved unconditionally at the end of Functions/CNMFE.m, and
% still carries neuron.Cn / the ring-model background either way).

Repo_Path = 'C:\Users\afermigier\Documents\GitHub\1Photon_Analysis';
addpath(genpath(Repo_Path));

mouse_color = '840R';
task        = 'EPM';

base_dir      = 'F:\Inscopix_Projects\202508_DualColorMiniscope';
channel_label = [mouse_color '_' task];
out_folder    = fullfile(base_dir, mouse_color, channel_label, 'processed_data');

if contains(mouse_color, 'R')
    mc_mat_path = fullfile(out_folder, 'MC.mat');
else
    mc_mat_path = fullfile(out_folder, 'MC_meanpad.mat');
end

%% Load the CNMFE result for this session
% Prefer the post-processed result if it exists (it's been through
% overlap-based cleanup); fall back to the raw CNMFE output otherwise -
% this is the file that exists even when 0 neurons were auto-detected.
postproc_path = fullfile(out_folder, 'Result_CNMFE_postprocessed.mat');
raw_path      = fullfile(out_folder, 'Result_CNMFE.mat');

if isfile(postproc_path)
    load(postproc_path, 'neuron');
else
    load(raw_path, 'neuron');
end

fprintf('Loaded neuron object: %d existing neuron(s)\n', size(neuron.A, 2));

%% Draw missed neurons on CNMFE's own background image
% neuron.Cn is the exact correlation image CNMFE itself used (see
% Outlines_Neurons.png) - use it as the drawing background so new ROIs
% line up with the same pixel grid CNMFE worked on. If you also want the
% PNR-weighted version CNMFE displays by default, use neuron.Cn .*
% max(neuron.PNR, 0) instead - both are stored on the object.
bg_img = neuron.Cn .* max(neuron.PNR, 0);
clim = prctile(bg_img(:), [1 99.5]);

ROIs = draw_polygon_rois(bg_img, clim);

%% Extract traces for the new ROIs the same way CNMFE extracts everyone else's
neuron = add_manual_neurons_cnmfe(neuron, mc_mat_path, ROIs);

%% Save the updated neuron object, and re-show the contours including the new ones
if isfile(postproc_path)
    save(postproc_path, 'neuron');
else
    save(raw_path, 'neuron');
end

Coor = neuron.show_contours(0.6);
saveas(gcf, fullfile(out_folder, 'Outlines_Neurons_with_manual.png'));

fprintf('Done - %d total neuron(s) after manual additions.\n', size(neuron.A, 2));

%% Optional: re-export to the same CNMFE_Data.mat / CSV format used downstream
Data.A = full(neuron.A);
Data.C = neuron.C;
Data.C_raw = neuron.C_raw;
Data.S = neuron.S;
Data.ids = neuron.ids;
save(fullfile(out_folder, 'CNMFE_Data_with_manual.mat'), 'Data');

out_csv = fullfile(out_folder, [channel_label '_Traces_with_manual.csv']);
T = size(neuron.C_raw, 2);
K = size(neuron.C_raw, 1);
var_names = arrayfun(@(k) sprintf('Neuron_%d', neuron.ids(k)), 1:K, 'UniformOutput', false);
Tbl = array2table(neuron.C_raw', 'VariableNames', var_names);
Tbl.Frame = (1:T)';
Tbl = Tbl(:, ['Frame', var_names]);
writetable(Tbl, out_csv);
fprintf('Saved traces (C_raw, all neurons) to %s\n', out_csv);