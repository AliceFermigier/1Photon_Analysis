function [S1, S2, Sxy_h, Sxy_v] = accumulate_local_corr_sums(frames, Bg_Sigma)
%ACCUMULATE_LOCAL_CORR_SUMS  Apply per-frame local-background subtraction
%   and accumulate the running sums needed for a local pixel-correlation
%   image, over an in-memory block of frames.
%
%   [S1, S2, Sxy_h, Sxy_v] = accumulate_local_corr_sums(frames, Bg_Sigma)
%
%   frames   - [H x W x n] block of frames (numeric, will be cast to double)
%   Bg_Sigma - sigma (pixels) of the Gaussian used to estimate each
%              frame's local background, subtracted before accumulating
%
%   See compute_correlation_image.m for what Bg_Sigma controls and why.

    [H, W, n] = size(frames);

    hsize = 2*ceil(3*Bg_Sigma) + 1;
    bg_kernel = fspecial('gaussian', hsize, Bg_Sigma);

    S1    = zeros(H, W);
    S2    = zeros(H, W);
    Sxy_h = zeros(H, W-1);
    Sxy_v = zeros(H-1, W);

    for f = 1:n
        frame = double(frames(:, :, f));
        bg = imfilter(frame, bg_kernel, 'replicate');
        Xf = frame - bg;   % spatially high-pass filtered frame

        S1 = S1 + Xf;
        S2 = S2 + Xf.^2;
        Sxy_h = Sxy_h + Xf(:, 1:end-1) .* Xf(:, 2:end);
        Sxy_v = Sxy_v + Xf(1:end-1, :) .* Xf(2:end, :);
    end
end