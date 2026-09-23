function [MaxProj_raw, MaxProj_corrected, GlobalIntensity, Correction] = ...
    compute_max_projection(mc_mat_path, Chunk_size, Crop_Margin, Correction_Smoothing_Window)
%COMPUTE_MAX_PROJECTION  Max-intensity projection over an entire MC.mat
%   video, with an optional global-intensity correction pass to remove
%   whole-FOV brightness swings (e.g. ambient light contamination from
%   the mouse moving relative to a light source) before taking the max,
%   so those swings don't dominate/wash out real cell contrast.
%
%   [MaxProj_raw, MaxProj_corrected, GlobalIntensity, Correction] = ...
%       compute_max_projection(mc_mat_path, Chunk_size, Crop_Margin, Correction_Smoothing_Window)
%
%   mc_mat_path   - path to a motion-corrected MC.mat (variable 'MC', [H W N])
%   Chunk_size    - frames processed at once, e.g. 3000
%   Crop_Margin   - pixels to exclude from each edge when computing the
%                   per-frame global intensity (optional, default 0).
%                   Use this to exclude the motion-correction border
%                   region - a reasonable value is
%                   ceil(max(abs(Shifts(:)))) from *_Shifts_applied.mat.
%                   Does NOT crop the output projections, only the region
%                   used to estimate global brightness.
%   Correction_Smoothing_Window - frames (optional, default 1 = no
%                   smoothing, i.e. every frame is normalized fully to
%                   the same global level). Raise this (e.g. to a few
%                   hundred frames) if you want to preserve slower,
%                   possibly-real, whole-FOV brightness trends and only
%                   remove faster ambient-light flicker. For building an
%                   anatomical max projection (this function's purpose),
%                   the default of full per-frame correction is usually
%                   what you want.
%
%   Outputs:
%   MaxProj_raw        - [H x W] max projection of the untouched data
%   MaxProj_corrected  - [H x W] max projection after per-frame global
%                         intensity correction
%   GlobalIntensity    - [N x 1] per-frame median intensity (interior
%                         region only, if Crop_Margin > 0)
%   Correction         - [N x 1] additive correction applied per frame

    if nargin < 3 || isempty(Crop_Margin)
        Crop_Margin = 0;
    end
    if nargin < 4 || isempty(Correction_Smoothing_Window)
        Correction_Smoothing_Window = 1;
    end

    MC_f = matfile(mc_mat_path);
    sz = size(MC_f, 'MC');
    H = sz(1); W = sz(2); N = sz(3);

    rows = (Crop_Margin+1):(H-Crop_Margin);
    cols = (Crop_Margin+1):(W-Crop_Margin);

    starts = 1:Chunk_size:N;

    % --- Pass 1: per-frame global intensity + raw max projection ---
    % (folded into one pass since both need to read every chunk anyway)
    GlobalIntensity = zeros(N, 1);
    MaxProj_raw = zeros(H, W);

    for c = 1:numel(starts)
        s0 = starts(c);
        n  = min(Chunk_size, N - s0 + 1);
        chunk = double(MC_f.MC(:, :, s0:s0+n-1));

        interior = chunk(rows, cols, :);
        GlobalIntensity(s0:s0+n-1) = median(reshape(interior, [], n), 1)';

        MaxProj_raw = max(MaxProj_raw, max(chunk, [], 3));
        fprintf('  [pass 1/2] frames %d-%d / %d\n', s0, s0+n-1, N);
    end

    % --- Compute per-frame correction ---
    if Correction_Smoothing_Window > 1
        Trend = movmedian(GlobalIntensity, Correction_Smoothing_Window);
    else
        Trend = GlobalIntensity;   % no smoothing: full per-frame normalization
    end
    target_level = median(Trend);
    Correction = target_level - Trend;   % additive, per frame

    % --- Pass 2: apply correction, accumulate corrected max projection ---
    MaxProj_corrected = zeros(H, W);
    for c = 1:numel(starts)
        s0 = starts(c);
        n  = min(Chunk_size, N - s0 + 1);
        chunk = double(MC_f.MC(:, :, s0:s0+n-1));

        corr_chunk = chunk + reshape(Correction(s0:s0+n-1), 1, 1, n);
        corr_chunk = max(0, min(65535, corr_chunk));

        MaxProj_corrected = max(MaxProj_corrected, max(corr_chunk, [], 3));
        fprintf('  [pass 2/2] frames %d-%d / %d\n', s0, s0+n-1, N);
    end
end