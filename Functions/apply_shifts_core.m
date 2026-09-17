function apply_shifts_core(raw_tif, Shifts, out_folder, out_filename, Chunk_size, Fill_Value)
%APPLY_SHIFTS_CORE  Apply an already frame-matched shift array to a raw tif stack.
%
%   apply_shifts_core(raw_tif, Shifts, out_folder, out_filename, Chunk_size, Fill_Value)
%
%   raw_tif      - path to the raw .tif to be corrected
%   Shifts       - [N x 2] array of [row_shift col_shift], one row per
%                  frame of raw_tif, already matched/ordered (no
%                  interleave logic happens here - that's the caller's job)
%   out_folder   - folder to write out_filename and Shifts_applied.mat into
%   out_filename - e.g. 'MC.mat' or 'MC_meanpad.mat'
%   Chunk_size   - frames processed at once, e.g. 3000
%   Fill_Value   - 'frame_mean' (default) | 'zero' | numeric constant
%                  see apply_shifts_to_red.m for details
%
%   Requires imtranslate_old.m, loadtiff.m, savefast.m from
%   fmi-basel/1Photon_Analysis on the MATLAB path.

    if nargin < 6 || isempty(Fill_Value)
        Fill_Value = 'frame_mean';
    end

    info = imfinfo(raw_tif);
    N = numel(info);
    H = info(1).Height;
    W = info(1).Width;

    if size(Shifts, 1) ~= N
        error('Shifts has %d rows but %s has %d frames - these must match exactly here.', ...
              size(Shifts,1), raw_tif, N);
    end

    if ~isfolder(out_folder)
        mkdir(out_folder);
    end
    out_mat = fullfile(out_folder, out_filename);

    MC = zeros(H, W, N, 'uint16'); %#ok<NASGU>
    savefast(out_mat, 'MC');
    clear MC
    MC_f = matfile(out_mat, 'Writable', true);

    [~, name_noext] = fileparts(out_filename);
    savefast(fullfile(out_folder, [name_noext '_Shifts_applied.mat']), 'Shifts');

    starts = 1:Chunk_size:N;
    for c = 1:numel(starts)
        s0 = starts(c);
        n  = min(Chunk_size, N - s0 + 1);

        chunk = loadtiff(raw_tif, s0, n);          % [H x W x n], raw uint16
        chunk_shifted = zeros(H, W, n, 'uint16');

        for f = 1:n
            frame_shift = Shifts(s0 + f - 1, :);   % [row_shift col_shift]
            frame_d = double(chunk(:,:,f));

            if ischar(Fill_Value) || isstring(Fill_Value)
                switch lower(char(Fill_Value))
                    case 'frame_mean'
                        F = mean(frame_d(:));
                    case 'zero'
                        F = 0;
                    otherwise
                        error('Unknown Fill_Value option: %s', Fill_Value);
                end
            else
                F = Fill_Value;   % fixed numeric constant
            end

            shifted = imtranslate_old(frame_d, frame_shift, F);
            chunk_shifted(:,:,f) = uint16(max(0, min(65535, round(shifted))));
        end

        MC_f.MC(:,:, s0:s0+n-1) = chunk_shifted;
        fprintf('  frames %d-%d / %d done\n', s0, s0+n-1, N);
    end

    fprintf('Saved motion-corrected channel to: %s\n', out_mat);
end