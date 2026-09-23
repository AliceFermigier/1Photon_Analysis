function Cn = compute_correlation_image(mc_mat_path, Chunk_size, Bg_Sigma)
%COMPUTE_CORRELATION_IMAGE  Local pixel-correlation image over an entire
%   MC.mat video - the same underlying concept CNMFE uses as its
%   background for show_contours (Cn), computed independently here so it
%   can be used as the background for manual polygon ROI drawing.
%
%   Cn = compute_correlation_image(mc_mat_path, Chunk_size, Bg_Sigma)
%
%   mc_mat_path - path to a motion-corrected MC.mat (variable 'MC', [H W N])
%   Chunk_size  - frames processed at once, e.g. 3000
%   Bg_Sigma    - sigma (pixels) of the Gaussian used to estimate each
%                 frame's local background, which is subtracted before
%                 computing neighbor correlations (optional, default 15).
%                 Set this a bit larger than a neuron's radius: it needs
%                 to be small enough to preserve cell-sized structure but
%                 large enough to remove broad, whole-FOV brightness
%                 patterns (ambient light contamination, vignetting).
%
%   Returns Cn, an [H x W] map where each pixel's value is the average
%   Pearson correlation, over time, between that pixel's (locally
%   background-subtracted) intensity and its 4-connected neighbors'.
%   Real cells produce spatially correlated fluctuations across their
%   whole soma (shared calcium dynamics), so they appear bright here;
%   independent pixel noise and background do not, so they stay low -
%   regardless of absolute brightness.
%
%   WHY THIS IS ROBUST TO GLOBAL ILLUMINATION DRIFT WITHOUT NEEDING A
%   SEPARATE TEMPORAL CORRECTION: subtracting each pixel's own local
%   Gaussian-blurred neighborhood, frame by frame, exactly cancels any
%   spatially broad/uniform brightness offset in that frame - no matter
%   how that offset varies from frame to frame. This is a fundamentally
%   different (and more robust) mechanism than the additive per-frame
%   normalization used for a raw max-intensity projection, so no
%   additional temporal global-intensity correction is needed here.
%   (It assumes the illumination artifact is spatially broad relative to
%   Bg_Sigma - true for ambient/room lighting - not a small, cell-sized
%   hotspot, which this would not remove.)

    if nargin < 3 || isempty(Bg_Sigma)
        Bg_Sigma = 15;
    end

    MC_f = matfile(mc_mat_path);
    sz = size(MC_f, 'MC');
    H = sz(1); W = sz(2); N = sz(3);

    % Running sums for a streaming Pearson correlation, accumulated
    % chunk-wise (via accumulate_local_corr_sums.m) so the full movie
    % never needs to be held in memory at once, then converted to a
    % correlation map at the end (via sums_to_corr_image.m) - both shared
    % with preview_correlation_image.m so the math is always identical.
    S1    = zeros(H, W);
    S2    = zeros(H, W);
    Sxy_h = zeros(H, W-1);
    Sxy_v = zeros(H-1, W);

    starts = 1:Chunk_size:N;
    for c = 1:numel(starts)
        s0 = starts(c);
        n  = min(Chunk_size, N - s0 + 1);
        chunk = MC_f.MC(:, :, s0:s0+n-1);

        [s1, s2, sxy_h, sxy_v] = accumulate_local_corr_sums(chunk, Bg_Sigma);
        S1 = S1 + s1; S2 = S2 + s2; Sxy_h = Sxy_h + sxy_h; Sxy_v = Sxy_v + sxy_v;

        fprintf('  [correlation image] frames %d-%d / %d\n', s0, s0+n-1, N);
    end

    Cn = sums_to_corr_image(S1, S2, Sxy_h, Sxy_v, N);
end