%% Step 4b: Run CNMFE on ALL red-channel sessions, all mice
% Discovers every {mouse}R / {mouse}R_{task} session under base_dir that
% has a motion-corrected MC.mat, and runs Functions/CNMFE.m on each one
% IN ITS OWN try/catch - so one bad/corrupt session logs an error and
% the batch keeps going, instead of the whole run crashing partway
% through (Functions/CNMFE.m's own internal loop has no such
% protection, which is why sessions are run one at a time here rather
% than handing it the whole list at once).
%
% This can take a very long time across many sessions (CNMFE is slow
% per-session even on one machine) - run it somewhere it can be left
% running unattended (overnight/over a weekend), and check
% cnmfe_batch_log.csv afterwards (or partway through - it's written
% incrementally, one row per session, as each one finishes).

Repo_Path = 'C:\Users\afermigier\Documents\GitHub\1Photon_Analysis';
addpath(genpath(Repo_Path));

base_dir = 'F:\Inscopix_Projects\202508_DualColorMiniscope';
Fs = 20;             % recording frame rate (frames/sec)
NW2Test = [];        % [] = let Functions/CNMFE.m pick worker counts automatically
Overwrite = false;   % if true, re-run sessions that already have a Result_CNMFE.mat

log_path = fullfile(base_dir, 'cnmfe_batch_log.csv');

%% Discover all {mouse}R mouse folders
all_entries = dir(base_dir);
mouse_dirs = all_entries([all_entries.isdir] & endsWith({all_entries.name}, 'R'));

if isempty(mouse_dirs)
    error('No "*R" mouse folders found directly under %s.', base_dir);
end

fprintf('Found %d red-channel mouse folder(s): %s\n', numel(mouse_dirs), strjoin({mouse_dirs.name}, ', '));

%% Build the full list of sessions to process
sessions = {};   % each entry: struct('folder', ..., 'mouse_color', ..., 'task', ...)
for m = 1:numel(mouse_dirs)
    mouse_color = mouse_dirs(m).name;
    mouse_path = fullfile(base_dir, mouse_color);

    session_entries = dir(fullfile(mouse_path, [mouse_color '_*']));
    session_entries = session_entries([session_entries.isdir]);

    for s = 1:numel(session_entries)
        session_name = session_entries(s).name;   % e.g. '840R_EPM'
        task = erase(session_name, [mouse_color '_']);
        session_folder = fullfile(mouse_path, session_name);

        sessions{end+1} = struct( ...    %#ok<SAGROW>
            'folder', session_folder, ...
            'mouse_color', mouse_color, ...
            'task', task);
    end
end

fprintf('Found %d session(s) total across all red-channel mice.\n\n', numel(sessions));

%% Initialize the log file (append if it already exists from a previous partial run)
if ~isfile(log_path)
    fid = fopen(log_path, 'w');
    fprintf(fid, 'mouse_color,task,status,message,elapsed_min,n_neurons,timestamp\n');
    fclose(fid);
end

%% Process each session
for i = 1:numel(sessions)
    sess = sessions{i};
    mc_mat_path = fullfile(sess.folder, 'processed_data', 'MC.mat');
    result_path = fullfile(sess.folder, 'processed_data', 'Result_CNMFE.mat');

    fprintf('[%d/%d] %s_%s ... ', i, numel(sessions), sess.mouse_color, sess.task);

    if ~isfile(mc_mat_path)
        fprintf('SKIPPED (no MC.mat)\n');
        append_log(log_path, sess, 'skipped', 'no MC.mat found', NaN, NaN);
        continue;
    end

    if isfile(result_path) && ~Overwrite
        fprintf('SKIPPED (already has Result_CNMFE.mat)\n');
        append_log(log_path, sess, 'skipped', 'already processed', NaN, NaN);
        continue;
    end

    tic;
    try
        CNMFE({sess.folder}, Fs, NW2Test);
        elapsed_min = toc / 60;

        load(result_path, 'neuron');
        n_neurons = size(neuron.A, 2);

        outline_png = fullfile(sess.folder, 'processed_data', 'Outlines_Neurons.png');
        if n_neurons > 0
            % Only rely on show_contours when there's something to draw -
            % it's untested bundled code for the zero-neuron edge case.
            neuron.show_contours(0.6);
            title(sprintf('%s_%s - %d neurons', sess.mouse_color, sess.task, n_neurons), 'Interpreter', 'none');
        else
            figure('Name', 'No neurons found');
            imagesc(neuron.Cn, prctile(neuron.Cn(:), [1 99.5])); axis image off; colormap gray;
            title(sprintf('%s_%s - 0 neurons (background shown for reference)', sess.mouse_color, sess.task), 'Interpreter', 'none');
        end
        saveas(gcf, outline_png);
        close(gcf);

        fprintf('OK (%d neurons, %.1f min)\n', n_neurons, elapsed_min);
        append_log(log_path, sess, 'ok', '', elapsed_min, n_neurons);
    catch ME
        elapsed_min = toc / 60;
        fprintf('FAILED: %s\n', ME.message);
        append_log(log_path, sess, 'failed', ME.message, elapsed_min, NaN);
    end
end

fprintf('\nBatch complete. See %s for the full log.\n', log_path);

%% ------------------------------------------------------------------
function append_log(log_path, sess, status, message, elapsed_min, n_neurons)
    fid = fopen(log_path, 'a');
    % escape commas/quotes in the message so the CSV doesn't break
    message = strrep(message, '"', '''');
    message = strrep(message, ',', ';');
    fprintf(fid, '%s,%s,%s,"%s",%.2f,%d,%s\n', ...
        sess.mouse_color, sess.task, status, message, elapsed_min, n_neurons, ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fclose(fid);
end