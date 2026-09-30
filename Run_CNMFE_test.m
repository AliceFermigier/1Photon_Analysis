%% Step 4a: Run CNMFE on ONE test session before batching
% Run this first on a single red-channel session to sanity-check
% parameters (see the tuning guidance already discussed - min_corr,
% min_pnr, min_corr_res, min_pnr_res, and the residual re-seeding patch
% to Functions/CNMFE.m) before committing to the full batch run in
% Step4_Run_CNMFE_batch_red.m. A single session can take a long time
% (potentially hours, depending on length/resolution/machine) - this is
% normal for CNMFE, not a hang.

Repo_Path = 'C:\Users\jcourtin\Documents\GitHub\Alice\1Photon_Analysis';
addpath(genpath(Repo_Path));

base_dir    = 'E:\Inscopix_Projects\202508_DualColorMiniscope';
mouse_color = '840R';
task        = 'EPM';

Fs = 20;          % recording frame rate (frames/sec) - same value used
                  % elsewhere in this pipeline (Recording_Speed)
NW2Test = [];     % [] lets Functions/CNMFE.m pick worker counts to try
                  % automatically (see its own NW2Testa logic) - only
                  % override this if you already know what works on
                  % your machine, e.g. NW2Test = [4 2 0]

%% Locate the session and check it's ready for CNMFE
channel_label  = [mouse_color '_' task];
session_folder = fullfile(base_dir, mouse_color, channel_label);
mc_mat_path    = fullfile(session_folder, 'processed_data', 'MC.mat');

if ~isfile(mc_mat_path)
    error('MC.mat not found at %s - run Step 1 (motion correction) for this session first.', mc_mat_path);
end

fprintf('Running CNMFE on: %s\n', session_folder);
fprintf('(Fs = %d Hz)\n', Fs);

%% Run CNMFE (Functions/CNMFE.m takes a CELL ARRAY of session folders)
tic;
CNMFE({session_folder}, Fs, NW2Test);
elapsed = toc;
fprintf('\nCNMFE finished in %.1f minutes.\n', elapsed/60);

%% Quick check of what it found
result_path = fullfile(session_folder, 'processed_data', 'Result_CNMFE.mat');
load(result_path, 'neuron');
fprintf('Found %d neuron(s).\n', size(neuron.A, 2));

outline_png = fullfile(session_folder, 'processed_data', 'Outlines_Neurons.png');

if size(neuron.A, 2) > 0
    % show_contours/plot_contours SHOULD also handle the zero-neuron
    % case gracefully (its loops just run zero times), but that's an
    % untested edge case in bundled code we didn't write - only rely on
    % it when there's actually something to draw. Below, the
    % zero-neuron branch saves the plain background directly instead,
    % which has no such dependency and always succeeds.
    neuron.show_contours(0.6);
    title(sprintf('%s - %d neurons', channel_label, size(neuron.A, 2)), 'Interpreter', 'none');
else
    fprintf('No neurons passed threshold - this is expected on some red-channel sessions.\n');
    fprintf('You can still use Step3_Add_Missed_Neurons.m to draw ROIs manually on neuron.Cn.\n');

    figure('Name', 'No neurons found');
    imagesc(neuron.Cn, prctile(neuron.Cn(:), [1 99.5])); axis image off; colormap gray;
    title(sprintf('%s - 0 neurons (background shown for reference)', channel_label), 'Interpreter', 'none');
end

saveas(gcf, outline_png);
fprintf('Saved image to %s\n', outline_png);