function apply_shifts_to_red(green_tif, shifts_mat, red_tif, out_folder, Chunk_size, Interleave_Mode)
%APPLY_SHIFTS_TO_RED  Apply a green channel's MC_Shifts.mat to the raw red channel.
%
%   apply_shifts_to_red(green_tif, shifts_mat, red_tif, out_folder, Chunk_size, Interleave_Mode)
%
%   green_tif       - path to the raw green .tif (used only for a frame-count sanity check)
%   shifts_mat      - path to the green session's processed_data\MC_Shifts.mat
%   red_tif         - path to the raw red .tif to be corrected
%   out_folder      - folder to write MC.mat and MC_Shifts_applied.mat into
%                      (typically {redsession}\processed_data)
%   Chunk_size      - frames processed at once, e.g. 3000
%   Interleave_Mode - 'auto' | 'green_first' | 'red_first' | 'same'
%                        green_first: pattern G,R,G,R,...,G  (N_green = N_red + 1)
%                        red_first:   pattern R,G,R,G,...,R  (N_red = N_green + 1)
%                        same:        N_green == N_red, matched 1:1
%
%   Requires imtranslate_old.m, loadtiff.m, savefast.m from
%   fmi-basel/1Photon_Analysis on the MATLAB path.

    % --- Load and concatenate the green channel's per-frame shifts ---
    % Shift_collection is a cell array (one cell per processing chunk),
    % each cell an [n_frames x 2] matrix of [row_shift col_shift], already
    % in the raw frame order. Simple vertical concatenation restores the
    % full-length, frame-ordered shift array.
    S = load(shifts_mat, 'Shift_collection');
    Shifts_green = cell2mat(S.Shift_collection(:));   % [N_green x 2]
    N_green = size(Shifts_green, 1);

    % --- Frame counts / dimensions ---
    info_G = imfinfo(green_tif);
    N_green_raw = numel(info_G);
    if N_green_raw ~= N_green
        warning(['Green raw tif has %d frames but MC_Shifts.mat has %d rows. ' ...
                 'This usually means a non-default T_DS_factor was used, or the ' ...
                 'session was reprocessed since. Proceeding with the shift array ' ...
                 'as-is - double check this session''s parameters.'], ...
                 N_green_raw, N_green);
    end

    info_R = imfinfo(red_tif);
    N_red  = numel(info_R);
    H = info_R(1).Height;
    W = info_R(1).Width;

    % --- Resolve the interleave offset between the two channels ---
    offset = N_green - N_red;
    switch Interleave_Mode
        case 'auto'
            if offset == 1
                mode = 'green_first';
            elseif offset == -1
                mode = 'red_first';
            elseif offset == 0
                mode = 'same';
            else
                error(['Unexpected frame count mismatch: %d green vs %d red frames ' ...
                       '(difference of %d, expected 0 or 1). Check acquisition/' ...
                       'deinterleaving for this session before proceeding.'], ...
                       N_green, N_red, offset);
            end
        otherwise
            mode = Interleave_Mode;
    end

    switch mode
        case 'green_first'   % pattern: G R G R ... G  (N_green = N_red + 1)
            Shifts_red = Shifts_green(1:N_red, :);
        case 'red_first'     % pattern: R G R G ... R  (N_red = N_green + 1)
            Shifts_red = Shifts_green(1:min(N_green, N_red), :);
            if N_red > size(Shifts_red, 1)
                warning('Extra unmatched red frame(s) at the end; reusing the final available shift for them.');
                Shifts_red = [Shifts_red; repmat(Shifts_green(end,:), N_red - size(Shifts_red,1), 1)];
            end
        case 'same'
            Shifts_red = Shifts_green(1:N_red, :);
        otherwise
            error('Unknown Interleave_Mode: %s', mode);
    end
    fprintf('Matched %d red frames to green shifts using mode: %s\n', N_red, mode);

    % --- Prepare output (.mat, same convention as the pipeline's own MC.mat) ---
    if ~isfolder(out_folder)
        mkdir(out_folder);
    end
    out_mat = fullfile(out_folder, 'MC.mat');

    MC = zeros(H, W, N_red, 'uint16'); %#ok<NASGU>
    savefast(out_mat, 'MC');
    clear MC
    MC_f = matfile(out_mat, 'Writable', true);

    savefast(fullfile(out_folder, 'MC_Shifts_applied.mat'), 'Shifts_red');

    % --- Apply shifts to the raw red channel, chunk by chunk ---
    starts = 1:Chunk_size:N_red;
    for c = 1:numel(starts)
        s0 = starts(c);
        n  = min(Chunk_size, N_red - s0 + 1);

        chunk = loadtiff(red_tif, s0, n);          % [H x W x n], raw uint16
        chunk_shifted = zeros(H, W, n, 'uint16');

        for f = 1:n
            frame_shift = Shifts_red(s0 + f - 1, :);   % [row_shift col_shift]
            % Same subpixel translation call the pipeline itself uses,
            % applied directly to the native-resolution raw frame.
            shifted = imtranslate_old(double(chunk(:,:,f)), frame_shift);
            chunk_shifted(:,:,f) = uint16(max(0, min(65535, round(shifted))));
        end

        MC_f.MC(:,:, s0:s0+n-1) = chunk_shifted;
        fprintf('  frames %d-%d / %d done\n', s0, s0+n-1, N_red);
    end

    fprintf('Saved motion-corrected red channel to: %s\n', out_mat);
end