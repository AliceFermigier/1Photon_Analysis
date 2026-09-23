function Traces = extract_roi_traces(mc_mat_path, ROIs, Chunk_size, out_csv_path)
%EXTRACT_ROI_TRACES  Extract per-frame mean fluorescence for each ROI from
%   a full MC.mat video and save to CSV. Uses the RAW (uncorrected)
%   motion-corrected data - the global intensity correction used for the
%   max projection is NOT applied here, since these traces are meant for
%   quantitative downstream analysis.
%
%   Traces = extract_roi_traces(mc_mat_path, ROIs, Chunk_size, out_csv_path)
%
%   mc_mat_path  - path to MC.mat (variable 'MC', [H W N])
%   ROIs         - struct array from draw_polygon_rois.m (needs .ID, .Mask)
%   Chunk_size   - frames processed at once, e.g. 3000
%   out_csv_path - where to save the CSV (Frame column + one column per ROI)
%
%   Returns Traces, an [N x n_rois] double array (same order as ROIs).

    MC_f = matfile(mc_mat_path);
    sz = size(MC_f, 'MC');
    H = sz(1); W = sz(2); N = sz(3);
    n_rois = numel(ROIs);

    roi_idx = cell(n_rois, 1);
    for r = 1:n_rois
        if ~isequal(size(ROIs(r).Mask), [H W])
            error('ROI %d mask size [%d %d] does not match video frame size [%d %d].', ...
                  ROIs(r).ID, size(ROIs(r).Mask,1), size(ROIs(r).Mask,2), H, W);
        end
        roi_idx{r} = find(ROIs(r).Mask);   % linear indices into an H x W frame
    end

    Traces = zeros(N, n_rois);
    starts = 1:Chunk_size:N;
    for c = 1:numel(starts)
        s0 = starts(c);
        n  = min(Chunk_size, N - s0 + 1);
        chunk = double(MC_f.MC(:, :, s0:s0+n-1));
        chunk = reshape(chunk, H*W, n);   % column-major linear index matches find(mask)

        for r = 1:n_rois
            Traces(s0:s0+n-1, r) = mean(chunk(roi_idx{r}, :), 1);
        end
        fprintf('  frames %d-%d / %d extracted\n', s0, s0+n-1, N);
    end

    var_names = [{'Frame'}, arrayfun(@(x) sprintf('Neuron_%d', x), [ROIs.ID], 'UniformOutput', false)];
    T = array2table([(1:N)', Traces], 'VariableNames', var_names);
    writetable(T, out_csv_path);
    fprintf('Saved fluorescence traces to %s\n', out_csv_path);
end