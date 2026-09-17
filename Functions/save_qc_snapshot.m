function save_qc_snapshot(raw_path, mc_mat_path, out_png_path, qc_start_frame, qc_num_frames, title_prefix)
%SAVE_QC_SNAPSHOT  Save a 2x2 QC figure (raw/corrected mean and std
%   projections) over a short window of frames, without popping up a
%   window (safe for unattended batch runs).
%
%   save_qc_snapshot(raw_path, mc_mat_path, out_png_path, qc_start_frame, qc_num_frames, title_prefix)
%
%   raw_path       - path to the raw .tif
%   mc_mat_path    - path to the corresponding MC.mat (variable 'MC')
%   out_png_path   - where to save the PNG
%   qc_start_frame - first frame of the QC window (1-indexed). Pick this
%                    a few seconds into the recording rather than frame 1,
%                    to avoid LED-onset illumination artifacts at the very
%                    start of acquisition.
%   qc_num_frames  - number of frames in the QC window, e.g. 200
%   title_prefix   - short label used in the figure title (e.g. session name)

    info = imfinfo(raw_path);
    N = numel(info);

    n = min(qc_num_frames, N - qc_start_frame + 1);
    if n < 1
        warning(['QC start frame %d is beyond the %d available frames in %s. ' ...
                 'Falling back to frame 1.'], qc_start_frame, N, raw_path);
        qc_start_frame = 1;
        n = min(qc_num_frames, N);
    end

    raw_chunk  = loadtiff(raw_path, qc_start_frame, n);
    MC_f       = matfile(mc_mat_path);
    corr_chunk = MC_f.MC(:, :, qc_start_frame:qc_start_frame+n-1);

    mean_raw  = mean(double(raw_chunk), 3);
    mean_corr = mean(double(corr_chunk), 3);
    std_raw   = std(double(raw_chunk), 0, 3);
    std_corr  = std(double(corr_chunk), 0, 3);

    clim_mean = [min([mean_raw(:); mean_corr(:)]), max([mean_raw(:); mean_corr(:)])];
    clim_std  = [min(prctile(std_raw(:),1), prctile(std_corr(:),1)), ...
                 max(prctile(std_raw(:),99), prctile(std_corr(:),99))];

    fig = figure('Visible', 'off', 'Position', [100 100 900 800]);

    ax1 = subplot(2,2,1); imagesc(mean_raw, clim_mean); axis image off; colormap(ax1, 'gray'); colorbar;
    title('Raw - mean');
    ax2 = subplot(2,2,2); imagesc(mean_corr, clim_mean); axis image off; colormap(ax2, 'gray'); colorbar;
    title('Corrected - mean');
    ax3 = subplot(2,2,3); imagesc(std_raw, clim_std); axis image off; colormap(ax3, 'hot'); colorbar;
    title('Raw - std');
    ax4 = subplot(2,2,4); imagesc(std_corr, clim_std); axis image off; colormap(ax4, 'hot'); colorbar;
    title('Corrected - std');

    % Overall title, drawn as a borderless textbox spanning the top of the
    % figure. Avoids sgtitle (MATLAB R2018b+ only, and can be unavailable
    % in some setups) so this works on any MATLAB version with no
    % additional dependencies.
    overall_title = sprintf('%s QC: frames %d-%d', title_prefix, qc_start_frame, qc_start_frame+n-1);
    annotation(fig, 'textbox', [0 0.955 1 0.04], ...
        'String', overall_title, 'Interpreter', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'FontSize', 11, 'EdgeColor', 'none');

    % exportgraphics is R2020a+; fall back to the older print() command if
    % it's unavailable, so this still works on earlier MATLAB versions.
    if exist('exportgraphics', 'file')
        exportgraphics(fig, out_png_path, 'Resolution', 150);
    else
        print(fig, out_png_path, '-dpng', '-r150');
    end
    close(fig);
    fprintf('Saved QC snapshot: %s\n', out_png_path);
end