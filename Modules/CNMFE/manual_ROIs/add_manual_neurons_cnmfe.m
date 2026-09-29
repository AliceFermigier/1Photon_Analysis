function neuron = add_manual_neurons_cnmfe(neuron, mc_mat_path, ROIs)
%ADD_MANUAL_NEURONS_CNMFE  Add manually-drawn polygon ROIs to a CNMFE
%   Sources2D result, extracting each one's trace THE SAME WAY CNMFE
%   extracts every other neuron's trace - by replicating the exact
%   per-neuron recipe in updateTemporal_endoscope.m:
%     1. subtract the estimated ring-model background AND the
%        reconstructed contribution of every already-existing neuron,
%        leaving a residual movie;
%     2. average the residual over the ROI's pixels (this is exactly
%        what CNMFE's own HALS update reduces to for an isolated,
%        non-overlapping spatial component - see the "temp = C(k,:) +
%        (U(k,:)-V(k,:)*C)/aa(k)" line in updateTemporal_endoscope.m);
%     3. estimate a baseline and noise level the same way
%        (estimate_baseline_noise.m + GetSn.m, keeping whichever gives
%        the lower noise estimate, exactly as CNMFE does);
%     4. deconvolve with deconvolveCa.m using the SAME deconv_options
%        the whole session was run with (neuron.options.deconv_options);
%     5. apply CNMFE's own A/C/C_raw/S normalization convention (A
%        scaled by the noise level, C/C_raw/S divided by it) so the new
%        neuron is stored in the exact same units as every other one.
%
%   neuron = add_manual_neurons_cnmfe(neuron, mc_mat_path, ROIs)
%
%   neuron      - a loaded Sources2D object (from Result_CNMFE.mat or
%                 Result_CNMFE_postprocessed.mat) - CAN have zero
%                 existing neurons (neuron.A empty), that's fine.
%   mc_mat_path - path to the SAME MC.mat CNMFE was actually run on
%                 (variable 'MC', [H W N])
%   ROIs        - struct array from draw_polygon_rois.m (fields ID,
%                 Vertices, Mask, Centroid) - draw these on neuron.Cn
%                 (or neuron.Cn .* max(neuron.PNR, 0), if you saved
%                 that combination) so they line up with CNMFE's own
%                 pixel grid.
%
%   Returns the same neuron object with each ROI appended as one more
%   column of A / row of C, C_raw, S, and one more entry in obj.ids.
%   ROIs are added one at a time, in order, so a later ROI's residual
%   already accounts for an earlier one added in the same call (in case
%   you draw two overlapping missed neurons together).
%
%   IMPORTANT CAVEATS:
%   - This treats each new ROI as spatially non-overlapping with every
%     EXISTING (already-detected) neuron for the purpose of the
%     averaging step in (2) - correct as long as CNMFE genuinely missed
%     that whole region (the usual reason you're adding it manually). If
%     a manually-added ROI heavily overlaps an existing neuron's
%     footprint, this simplified single-component update will not
%     properly deblend the two - in that case it's more correct to
%     rerun neuron.update_spatial_parallel / update_temporal_parallel on
%     the whole object after appending A, letting CNMFE's own joint fit
%     handle the overlap (slower, but exactly what the full pipeline
%     would do).
%   - neuron.reconstruct_background() needs neuron.P.mat_data to still
%     point to a valid, reachable memory-mapped data file from the
%     original run. If that's been moved/deleted, this falls back to
%     neuron.reconstruct_b0() (a constant per-pixel baseline only, no
%     ring-model background dynamics) and prints a warning - traces will
%     be less clean in that case.

    d1 = neuron.options.d1;
    d2 = neuron.options.d2;

    MC_f = matfile(mc_mat_path);
    sz = size(MC_f, 'MC');
    T = sz(3);
    fprintf('[add_manual_neurons_cnmfe] loading full video (%d x %d x %d)...\n', d1, d2, T);
    Y = double(MC_f.MC(:, :, :));
    Y = reshape(Y, d1*d2, T);

    fprintf('[add_manual_neurons_cnmfe] reconstructing CNMFE''s own background...\n');
    try
        Ybg = neuron.reconstruct_background([1, T]);
        Ybg = reshape(Ybg, d1*d2, T);
    catch ME
        warning(['reconstruct_background() failed (%s). Falling back to a constant ' ...
                 'per-pixel baseline (reconstruct_b0) - no ring-model dynamics. ' ...
                 'Traces for the manually added neuron(s) will be noisier than usual.'], ME.message);
        b0_ = neuron.reconstruct_b0();
        Ybg = repmat(reshape(b0_, d1*d2, 1), 1, T);
    end

    deconv_options_0 = neuron.options.deconv_options;

    for r = 1:numel(ROIs)
        mask = ROIs(r).Mask;
        idx = find(mask(:));
        aa = numel(idx);
        if aa == 0
            warning('ROI %d has an empty mask - skipping.', ROIs(r).ID);
            continue;
        end

        % Residual = raw - background - every EXISTING neuron's own
        % reconstructed contribution (including any manual neuron
        % already appended earlier in this same loop).
        Y_roi = Y(idx, :);
        Ybg_roi = Ybg(idx, :);
        if size(neuron.A, 2) > 0
            existing_contrib = full(neuron.A(idx, :)) * neuron.C;
        else
            existing_contrib = 0;
        end
        resid = Y_roi - Ybg_roi - existing_contrib;

        % Single-component HALS update, reduces to a plain mean over the
        % ROI's pixels of the residual (see updateTemporal_endoscope.m).
        temp = sum(resid, 1) / aa;

        % Baseline + noise estimate, exactly as updateTemporal_endoscope.m
        [b_hist, sn_hist] = estimate_baseline_noise(temp);
        b = mean(temp(temp < median(temp)));
        sn_psd = GetSn(temp);
        if sn_psd < sn_hist
            sn = sn_psd;
        else
            sn = sn_hist;
            b = b_hist;
        end
        temp = temp - b;

        % Deconvolution, using the SAME options the whole session used.
        if neuron.options.deconv_flag
            [ck, sk, deconv_options] = deconvolveCa(temp, deconv_options_0, 'maxIter', 2, 'sn', sn);
            c_raw_final = temp - deconv_options.b;
        else
            ck = max(0, temp);
            sk = zeros(size(temp));
            c_raw_final = temp;
        end

        % Append as one more component, using CNMFE's own normalization
        % convention (A carries the amplitude scale sn; C/C_raw/S are
        % stored divided by it).
        new_A_col = sparse(idx, ones(aa,1), sn, d1*d2, 1);
        new_C_row = ck / sn;
        new_Craw_row = c_raw_final / sn;
        new_S_row = sk / sn;

        neuron.A = [neuron.A, new_A_col];
        neuron.C = [neuron.C; new_C_row];
        neuron.C_raw = [neuron.C_raw; new_Craw_row];
        neuron.S = [neuron.S; new_S_row];
        if isempty(neuron.ids)
            new_id = 1;
        else
            new_id = max(neuron.ids) + 1;
        end
        neuron.ids(end+1) = new_id;

        fprintf('[add_manual_neurons_cnmfe] added manual neuron (new id %d, %d pixels, sn=%.3g)\n', ...
                new_id, aa, sn);
    end
end