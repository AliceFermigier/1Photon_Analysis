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
%       {Root}\{mouse}G\{mouse}G_{experiment}\{mouse}G_{experiment}.tiff
%       {Root}\{mouse}R\{mouse}R_{experiment}\{mouse}R_{experiment}.tiff
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
% Written by Alice Fermigier, building on fmi-basel/1Photon_Analysis
% (Julian Hinz, Luthi lab)

%% ============ USER SETTINGS ============

Repo_Path = 'C:\Users\jcourtin\Documents\GitHub\Alice\1Photon_Analysis';
addpath(genpath(Repo_Path));

Root_Folder = 'F:\Inscopix_Projects\DualColor_Deinterleaved_data';

Chunk_size = 3000;         % frames processed at once (matches pipeline default)
Interleave_Mode = 'auto';  % 'auto' | 'green_first' | 'red_first' | 'same'
                            %   green_first: pattern G,R,G,R,...,G  (N_green = N_red + 1)  <- your case
                            %   red_first:   pattern R,G,R,G,...,R  (N_red = N_green + 1)
                            %   same:        N_green == N_red, matched 1:1
Fill_Value = 'frame_mean';  % 'frame_mean' (recommended) | 'zero' | numeric constant
Overwrite = false;         % if false, skip sessions whose output already exists

Process_Green = true;       % also regenerate the green channel with the same Fill_Value,
                             % saved as MC_meanpad.mat (does NOT overwrite the pipeline's
                             % own zero-padded MC.mat)

Save_QC = true;              % save a quick QC snapshot (mean+std, raw vs corrected) per channel
Recording_Speed = 20;        % Hz, per-channel frame rate after deinterleaving (matches pipeline default)
QC_Start_Sec = 5;           % seconds into the recording where the QC window starts -
                              % kept away from frame 1 to avoid LED-onset illumination artifacts
QC_Num_Frames = 200;         % number of frames in the QC window
Overwrite_QC = false;        % if false, skip QC snapshots that already exist (independent of Overwrite)

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
 
        green_tif  = fullfile(green_session_path, [exp_folder_G '.tiff']);
        red_tif    = fullfile(red_session_path,   [exp_folder_R '.tiff']);
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
        else
            fprintf('\n=== Processing %s -> %s ===\n', exp_folder_G, exp_folder_R);
            apply_shifts_to_red(green_tif, shifts_mat, red_tif, out_folder, Chunk_size, Interleave_Mode, Fill_Value); %#ok<*NODEF> - function file on path
        end
 
        green_out_folder = fullfile(green_session_path, 'processed_data');
        if Process_Green
            green_meanpad_mat = fullfile(green_out_folder, 'MC_meanpad.mat');
            if isfile(green_meanpad_mat) && ~Overwrite
                fprintf('Already processed (mean-padded green): %s (skipping)\n', exp_folder_G);
            else
                fprintf('=== Regenerating green with Fill_Value = %s: %s ===\n', Fill_Value, exp_folder_G);
                apply_shifts_to_green(green_tif, shifts_mat, green_out_folder, Chunk_size, Fill_Value);
            end
        end
 
        % --- QC snapshots (mean+std, raw vs corrected), a bit into the recording ---
        if Save_QC
            QC_Skip_Frames = round(QC_Start_Sec * Recording_Speed);
 
            % Red channel QC
            red_qc_png = fullfile(out_folder, [exp_folder_R '_QC.png']);
            if isfile(red_qc_png) && ~Overwrite_QC
                fprintf('QC snapshot already exists: %s (skipping)\n', red_qc_png);
            elseif isfile(out_mat)
                save_qc_snapshot(red_tif, out_mat, red_qc_png, QC_Skip_Frames + 1, QC_Num_Frames, exp_folder_R);
            else
                warning('Skipping red QC snapshot for %s - MC.mat not found (processing may have failed or been skipped).', exp_folder_R);
            end
 
            % Green channel QC - only meaningful if the mean-padded green was (re)generated
            if Process_Green
                green_meanpad_mat = fullfile(green_out_folder, 'MC_meanpad.mat');
                green_qc_png = fullfile(green_out_folder, [exp_folder_G '_QC.png']);
                if isfile(green_qc_png) && ~Overwrite_QC
                    fprintf('QC snapshot already exists: %s (skipping)\n', green_qc_png);
                elseif isfile(green_meanpad_mat)
                    save_qc_snapshot(green_tif, green_meanpad_mat, green_qc_png, QC_Skip_Frames + 1, QC_Num_Frames, exp_folder_G);
                else
                    warning('Skipping green QC snapshot for %s - MC_meanpad.mat not found.', exp_folder_G);
                end
            end
        end
    end
end
 
fprintf('\nAll sessions done.\n');
 
% Core logic now lives in the standalone function file apply_shifts_to_red.m
% (must be on the MATLAB path, e.g. saved next to this script) so it can
% also be called directly for single-session testing - see
% Test_Apply_MC_to_Red_single_session.m