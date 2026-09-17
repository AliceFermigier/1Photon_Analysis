function apply_shifts_to_green(green_tif, shifts_mat, out_folder, Chunk_size, Fill_Value, out_filename)
%APPLY_SHIFTS_TO_GREEN  Regenerate the green channel's own motion correction
%   using apply_shifts_core, so it can use the same border Fill_Value as
%   the red channel (the original pipeline always used 'zero').
%
%   apply_shifts_to_green(green_tif, shifts_mat, out_folder, Chunk_size, Fill_Value, out_filename)
%
%   green_tif    - path to the raw green .tif
%   shifts_mat   - path to this session's processed_data\MC_Shifts.mat
%                  (the same file the original pipeline produced)
%   out_folder   - folder to write the output into (typically the same
%                  processed_data folder the pipeline already used)
%   Chunk_size   - frames processed at once, e.g. 3000
%   Fill_Value   - 'frame_mean' (default) | 'zero' | numeric constant
%   out_filename - output filename (optional, default 'MC_meanpad.mat').
%                  Deliberately NOT 'MC.mat' by default, so this doesn't
%                  silently overwrite the original pipeline's own
%                  zero-padded MC.mat - pass 'MC.mat' explicitly with
%                  Overwrite handling in the caller if you do want to
%                  replace it.
%
%   Requires imtranslate_old.m, loadtiff.m, savefast.m, apply_shifts_core.m
%   on the MATLAB path.

    if nargin < 5 || isempty(Fill_Value)
        Fill_Value = 'frame_mean';
    end
    if nargin < 6 || isempty(out_filename)
        out_filename = 'MC_meanpad.mat';
    end

    S = load(shifts_mat, 'Shift_collection');
    Shifts_green = cell2mat(S.Shift_collection(:));   % [N_green x 2]

    info_G = imfinfo(green_tif);
    N_green_raw = numel(info_G);
    if N_green_raw ~= size(Shifts_green, 1)
        warning(['Green raw tif has %d frames but MC_Shifts.mat has %d rows. ' ...
                 'This usually means a non-default T_DS_factor was used, or the ' ...
                 'session was reprocessed since. Proceeding with the shift array ' ...
                 'as-is - double check this session''s parameters.'], ...
                 N_green_raw, size(Shifts_green, 1));
    end

    apply_shifts_core(green_tif, Shifts_green, out_folder, out_filename, Chunk_size, Fill_Value);
end