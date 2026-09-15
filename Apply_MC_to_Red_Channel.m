%% Apply Green-Channel Motion Correction to the Red Channel
%
% Reuses the per-frame [row_shift col_shift] values that
% fmi-basel/1Photon_Analysis computed and saved in each green session's
% processed_data\MC_Shifts.mat, and applies the IDENTICAL translation to
% the corresponding raw red frame, using the exact same subpixel
% translation routine the pipeline itself uses (Utilities/imtranslate_old.m).
%
% Assumes:
%   - Interleaved acquisition (G,R,G,R,...) that was de-interleaved into
%     two separate tif stacks by the Inscopix software, so green and red
%     frame i were acquired ~half a frame period apart (close enough in
%     time that reusing green's shift for the matching red frame is a
%     good approximation - motion correction operates on a much slower
%     timescale than one inter-frame interval).
%   - Folder layout:
%       {Root}\{mouse}G\{mouse}G_{experiment}\{mouse}G_{experiment}.tif
%       {Root}\{mouse}R\{mouse}R_{experiment}\{mouse}R_{experiment}.tif
%   - The green session already has processed_data\MC_Shifts.mat
%     (i.e. you already ran the repo's oneP_Image_Analysis.m on the green
%     channel).
%
% Output:
%   {mouse}R_{experiment}\processed_data\MC.mat        - motion-corrected
%                                                         red channel, same
%                                                         convention as the
%                                                         pipeline's own MC.mat
%   {mouse}R_{experiment}\processed_data\MC_Shifts_applied.mat
%                                                       - the exact per-frame
%                                                         shifts used, for
%                                                         provenance
%
% Written by/for: [your name], building on fmi-basel/1Photon_Analysis
% (Julian Hinz, Luthi lab)

%% ============ USER SETTINGS ============

Repo_Path = 'C:\path\to\1Photon_Analysis';     % root of the cloned repo
addpath(genpath(Repo_Path));

Root_Folder = 'I:\Inscopix_Projects\DualColor_Deinterleaved_data';

Chunk_size = 3000;         % frames processed at once (matches pipeline default)
Interleave_Mode = 'auto';  % 'auto' | 'green_first' | 'red_first' | 'same'
                            %   green_first: pattern G,R,G,R,...,G  (N_green = N_red + 1)  <- your case
                            %   red_first:   pattern R,G,R,G,...,R  (N_red = N_green + 1)
                            %   same:        N_green == N_red, matched 1:1
Overwrite = false;         % if false, skip sessions whose MC.mat already exists

%% ============ FIND AND PROCESS SESSIONS ============

green_mouse_folders = dir(fullfile(Root_Folder, '*G'));
green_mouse_folders = green_mouse_folders([green_mouse_folders.isdir]);

for m = 1:numel(green_mouse_folders)

    mouse_G_name = green_mouse_folders(m).name;        % e.g. '839G'
    mouse_id     = mouse_G_name(1:end-1);               % e.g. '839'
    mouse_R_name = [mouse_id 'R'];
    mouse_G_path = fullfile(Root_Folder, mouse_G_name);
    mouse_R_path = fullfile(Root_Folder, mouse_R_name);

    if ~isfolder(mouse_R_path)
        warning('No matching red folder for %s, skipping.', mouse_G_name);
        continue
    end

    session_folders = dir(fullfile(mouse_G_path, [mouse_G_name '_*']));
    session_folders = session_folders([session_folders.isdir]);

    for s = 1:numel(session_folders)

        exp_folder_G = session_folders(s).name;                       % '839G_EPM'
        exp_name     = extractAfter(exp_folder_G, [mouse_G_name '_']); % 'EPM'
        exp_folder_R = [mouse_R_name '_' exp_name];                    % '839R_EPM'

        green_session_path = fullfile(mouse_G_path, exp_folder_G);
        red_session_path   = fullfile(mouse_R_path, exp_folder_R);

        green_tif  = fullfile(green_session_path, [exp_folder_G '.tif']);
        red_tif    = fullfile(red_session_path,   [exp_folder_R '.tif']);
        shifts_mat = fullfile(green_session_path, 'processed_data', 'MC_Shifts.mat');

        if ~isfolder(red_session_path) || ~isfile(red_tif)
            warning('Missing red session/tif for %s, skipping.', exp_folder_G);
            continue
        end
        if ~isfile(shifts_mat)
            warning('No MC_Shifts.mat for %s - run green motion correction first. Skipping.', exp_folder_G);
            continue
        end

        out_folder = fullfile(red_session_path, 'processed_data');
        out_mat    = fullfile(out_folder, 'MC.mat');
        if isfile(out_mat) && ~Overwrite
            fprintf('Already processed: %s (skipping, Overwrite = false)\n', exp_folder_R);
            continue
        end

        fprintf('\n=== Processing %s -> %s ===\n', exp_folder_G, exp_folder_R);
        apply_shifts_to_red(green_tif, shifts_mat, red_tif, out_folder, Chunk_size, Interleave_Mode); %#ok<*NODEF> - function file on path
    end
end

fprintf('\nAll sessions done.\n');