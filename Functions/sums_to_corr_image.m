function Cn = sums_to_corr_image(S1, S2, Sxy_h, Sxy_v, N)
%SUMS_TO_CORR_IMAGE  Convert accumulated per-pixel sums into a local
%   pixel-correlation image. Shared by compute_correlation_image.m (which
%   streams the sums chunk-wise over a full video) and
%   preview_correlation_image.m (which accumulates them once over a short
%   in-memory window, for fast parameter tuning) - kept as one function
%   so both always compute the correlation identically.
%
%   Cn = sums_to_corr_image(S1, S2, Sxy_h, Sxy_v, N)
%
%   S1    - [H x W] sum over frames of the (background-subtracted) pixel value
%   S2    - [H x W] sum over frames of that value squared
%   Sxy_h - [H x (W-1)] sum over frames of (pixel .* right-neighbor)
%   Sxy_v - [(H-1) x W] sum over frames of (pixel .* bottom-neighbor)
%   N     - number of frames the sums were accumulated over

    H = size(S1, 1);
    W = size(S1, 2);

    meanX = S1 / N;
    varX  = S2 / N - meanX.^2;
    varX(varX < 0) = 0;   % guard against tiny negative values from floating-point error

    % Horizontal edges (pixel <-> right neighbor)
    covXY_h = Sxy_h / N - meanX(:, 1:end-1) .* meanX(:, 2:end);
    denom_h = sqrt(varX(:, 1:end-1) .* varX(:, 2:end));
    corr_h = zeros(H, W-1);
    valid_h = denom_h > 0;
    corr_h(valid_h) = covXY_h(valid_h) ./ denom_h(valid_h);

    % Vertical edges (pixel <-> bottom neighbor)
    covXY_v = Sxy_v / N - meanX(1:end-1, :) .* meanX(2:end, :);
    denom_v = sqrt(varX(1:end-1, :) .* varX(2:end, :));
    corr_v = zeros(H-1, W);
    valid_v = denom_v > 0;
    corr_v(valid_v) = covXY_v(valid_v) ./ denom_v(valid_v);

    % Average the (up to 4) edge correlations incident on each pixel
    Cn_sum = zeros(H, W);
    Cn_cnt = zeros(H, W);

    Cn_sum(:, 1:end-1) = Cn_sum(:, 1:end-1) + corr_h;
    Cn_cnt(:, 1:end-1) = Cn_cnt(:, 1:end-1) + 1;
    Cn_sum(:, 2:end)   = Cn_sum(:, 2:end)   + corr_h;
    Cn_cnt(:, 2:end)   = Cn_cnt(:, 2:end)   + 1;

    Cn_sum(1:end-1, :) = Cn_sum(1:end-1, :) + corr_v;
    Cn_cnt(1:end-1, :) = Cn_cnt(1:end-1, :) + 1;
    Cn_sum(2:end, :)   = Cn_sum(2:end, :)   + corr_v;
    Cn_cnt(2:end, :)   = Cn_cnt(2:end, :)   + 1;

    Cn = Cn_sum ./ Cn_cnt;
end