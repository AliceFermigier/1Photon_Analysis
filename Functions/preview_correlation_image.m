function Cn_list = preview_correlation_image(mc_mat_path, Frame_Start, N_Frames, Bg_Sigma_List)
%PREVIEW_CORRELATION_IMAGE  Compute and display the local correlation
%   image for several candidate Bg_Sigma values, on a short in-memory
%   window of frames - for quickly tuning Bg_Sigma before committing to
%   compute_correlation_image.m on the full video.
%
%   Cn_list = preview_correlation_image(mc_mat_path, Frame_Start, N_Frames, Bg_Sigma_List)
%
%   mc_mat_path   - path to a motion-corrected MC.mat (variable 'MC', [H W N])
%   Frame_Start   - first frame of the preview window (1-indexed). Pick
%                   this a bit into the recording, same reasoning as the
%                   QC snapshots - avoids any LED-onset artifacts.
%   N_Frames      - number of frames to use, e.g. 500-1000. Doesn't need
%                   to be huge - just enough for a stable correlation
%                   estimate; this is purely for visual parameter tuning.
%   Bg_Sigma_List - vector of Bg_Sigma values to compare, e.g. [5 10 15 20 30]
%
%   Returns Cn_list, a cell array of the resulting [H x W] correlation
%   images, one per Bg_Sigma value, in the same order.
%
%   HOW TO USE THE RESULT: look for the smallest Bg_Sigma at which real
%   neurons appear as clean, well-separated, roughly round/soma-shaped
%   blobs of high correlation - not fragmented/dim (Bg_Sigma too small:
%   the local background estimate starts to include the cell's own
%   signal, partially subtracting it out) and not merged into large
%   blotchy patches or still showing broad illumination structure
%   (Bg_Sigma too large: it stops acting as a local background and
%   starts just smoothing/failing to remove broad brightness patterns).
%   A good starting point is roughly 1-1.5x the neuron radius in pixels -
%   see measure_neuron_radius.m to measure that directly on your data.

    MC_f = matfile(mc_mat_path);
    sz = size(MC_f, 'MC');
    N = sz(3);

    n = min(N_Frames, N - Frame_Start + 1);
    if n < N_Frames
        warning('Only %d frames available from Frame_Start = %d (requested %d).', n, Frame_Start, N_Frames);
    end

    frames = MC_f.MC(:, :, Frame_Start:Frame_Start+n-1);

    n_sigma = numel(Bg_Sigma_List);
    Cn_list = cell(n_sigma, 1);

    n_cols = ceil(sqrt(n_sigma));
    n_rows = ceil(n_sigma / n_cols);
    figure('Name', 'Bg_Sigma preview', 'Position', [100 100 300*n_cols 300*n_rows]);

    for k = 1:n_sigma
        Bg_Sigma = Bg_Sigma_List(k);
        [S1, S2, Sxy_h, Sxy_v] = accumulate_local_corr_sums(frames, Bg_Sigma);
        Cn = sums_to_corr_image(S1, S2, Sxy_h, Sxy_v, n);
        Cn_list{k} = Cn;

        subplot(n_rows, n_cols, k);
        imagesc(Cn, prctile(Cn(:), [1 99.5])); axis image off; colormap gray;
        title(sprintf('Bg\\_Sigma = %g', Bg_Sigma), 'Interpreter', 'tex');

        fprintf('Bg_Sigma = %g done (%d/%d)\n', Bg_Sigma, k, n_sigma);
    end
end