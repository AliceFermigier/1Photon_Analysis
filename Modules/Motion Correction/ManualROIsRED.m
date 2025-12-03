%% ManualROIsRED.m
% Manually extract cells from RED movie (expects processed_data/MC_red.mat)
% Saves:
%  - processed_data/RED_manual_results.mat      (polygons, masks, traces, dF/F, Red_MC)
%  - processed_data/RED_manual_CNMFE_style.mat  (A, C_raw, C, S, Cn, Df, C_df)
%  - processed_data/Outlines_Neurons_RED.png
%
% Usage: put All_nam, Data_Format_in, T_DS_factor, Spatial_Downsampling in workspace
% then run this script.

%% ---------- PARAMETERS (edit as desired) ----------
% Baseline/dF/F options
DF_PRCTILE = 10;       % percentile for baseline F0 when no window used (default 10)
DF_WINDOW = [];        % sliding window in frames (set e.g. 300 to use running percentile). [] = global

% Neuropil subtraction
DO_NEUROPIL = false;    % estimate neuropil and subtract (common for 1p)
NEUROPIL_RADIUS = 8;   % outer radius of neuropil ring (pixels)
NEUROPIL_INNER = 4;    % inner radius to exclude ROI area (pixels)
NEUROPIL_SCALE = 0.7;  % factor * neuropil_trace subtracted from ROI trace

% Visualization / filtering
gSig = 4;              % gaussian kernel sigma (same as Motion_correction)
psf = fspecial('gaussian', ceil(gSig*4 + 1), gSig);

% Save options
SAVE_OUTLINE_PNG = true;

%% ---------- Sanity checks ----------
if ~exist('All_nam','var')
    error('All_nam is not in workspace. Provide same All_nam used for motion correction.');
end
if ~exist('Data_Format_in','var') || ~exist('T_DS_factor','var') || ~exist('Spatial_Downsampling','var')
    error('Make sure Data_Format_in, T_DS_factor and Spatial_Downsampling are defined in the workspace.');
end

Resize_factor = Spatial_Downsampling; % keep naming consistent

%% ---------- Main loop over animals & sessions ----------
for p = 1:size(All_nam,1)
    subfolders = All_nam{p}; % cell array of session folder paths for this animal
    
    for pp = 1:numel(subfolders)
        Current_folder = [subfolders{pp} filesep];
        proc_dir = fullfile(Current_folder, 'processed_data');
        red_mat_path = fullfile(proc_dir, 'MC_red.mat');
        
        if ~exist(red_mat_path,'file')
            warning('No MC_red.mat found in %s — run ApplyShiftsRED first. Skipping.', Current_folder);
            continue
        end
        
        % try to create a sensible session name for figure titles and files
        [~, foldername] = fileparts(subfolders{pp});
        red_name = foldername;
        
        % ---------- Load corrected red movie ----------
        tmp = load(red_mat_path);
        if isfield(tmp,'Red_MC')
            Red_MC = tmp.Red_MC;   % H x W x T
        else
            % If user saved variable with different name (rare), try to detect 3D array
            fn = fieldnames(tmp);
            found = false;
            for k = 1:numel(fn)
                v = tmp.(fn{k});
                if ndims(v) == 3 && size(v,3) > 1
                    Red_MC = v;
                    found = true;
                    break
                end
            end
            if ~found
                error('MC_red.mat found but no 3D array inside to use as Red_MC.');
            end
        end
        
        [H, W, Tframes] = size(Red_MC);
        d = H*W;
        
        % ---------- Create filtered max projection for drawing ----------
        % Convert to double only once
        RedD = double(Red_MC);         % H x W x T  (may be memory heavy but necessary for filtering)
        
        % Apply gaussian filter frame-by-frame (imfilter accepts 3D directly with same kernel)
        % imfilter will apply filter spatially to each frame
        RedD_filt = imfilter(RedD, psf, 'replicate');
        
        % Max projection
        max_red = max(RedD_filt, [], 3);
        
        % Optional contrast enhance (choose as needed, currently commented)
        % max_red = imadjust(mat2gray(max_red));
        % max_red = adapthisteq(mat2gray(max_red));  % CLAHE - good for faint features
        
        hfig = figure('Name', ['Max RED Projection: ' red_name], 'NumberTitle','off');
        imagesc(max_red); axis image off; colormap gray;
        title(['Max RED Projection (filtered): ' red_name]);
        drawnow;
        
        % ---------- Manual ROI drawing ----------
        Red_ROIs = {};  % polygon vertex lists
        Red_Masks = {}; % logical masks
        keep_drawing = true;
        uiwait(msgbox('Draw ROIs on the RED max projection. Double-click to finish each ROI. Click OK to start.', 'Draw ROIs','modal'));
        
        while keep_drawing
            try
                h = drawpolygon('LineWidth',1.5,'Color','r');    % user draws polygon and double-clicks to close
            catch ME
                close(hfig);
                rethrow(ME);
            end
            pos = h.Position; % Nx2 [x y] (x: column, y: row)
            Red_ROIs{end+1} = pos;
            mask = poly2mask(pos(:,1), pos(:,2), H, W);
            Red_Masks{end+1} = mask;
            
            choice = questdlg('Add another ROI?', 'ROI','Yes','No','Yes');
            if strcmp(choice,'No'), keep_drawing = false; end
        end
        close(hfig);
        
        % ---------- Vectorized trace extraction ----------
        % Reshape movie to d x T (double) for matrix ops
        Yr = reshape(RedD, d, Tframes);   % large but needed for fast extraction
        clear RedD RedD_filt;  % free up memory if possible (max_red kept)
        
        numROIs = numel(Red_Masks);
        if numROIs == 0
            warning('No ROIs drawn for session %s. Skipping.', Current_folder);
            continue
        end
        
        A_sparse = sparse(d, numROIs);
        C_raw = zeros(numROIs, Tframes);
        ROI_pixel_counts = zeros(numROIs,1);
        
        for r = 1:numROIs
            m = Red_Masks{r};
            vec = reshape(m, [], 1);
            ROI_pixel_counts(r) = nnz(vec);
            if ROI_pixel_counts(r) == 0
                warning('ROI %d has zero pixels (skipped).', r);
                continue
            end
            A_sparse(:, r) = sparse(double(vec));   % d x 1 column
            % mean across pixels: (vec' * Yr) / sum(vec)
            C_raw(r, :) = (vec' * Yr) / ROI_pixel_counts(r);
        end
        
        % ---------- Optional neuropil subtraction ----------
        if DO_NEUROPIL
            fprintf('Estimating neuropil and subtracting (radius %d inner %d, scale %.2f)\n', NEUROPIL_RADIUS, NEUROPIL_INNER, NEUROPIL_SCALE);
            % Precompute structuring elements
            se_outer = strel('disk', NEUROPIL_RADIUS, 0);
            se_inner = strel('disk', NEUROPIL_INNER, 0);
            C_neuropil = zeros(numROIs, Tframes);
            
            % For morphological ops we need 2D masks
            for r = 1:numROIs
                roi_mask = Red_Masks{r};
                % dilate ROI to neuropil ring
                dil = imdilate(roi_mask, se_outer);
                inner = imdilate(roi_mask, se_inner);
                ring = dil & ~inner;
                % remove any overlap with other ROIs (exclude other ROI pixels)
                for rr = 1:numROIs
                    if rr == r, continue; end
                    ring(Red_Masks{rr}) = false;
                end
                if nnz(ring) < 10
                    % fallback if neuropil ring too small
                    ring = dil & ~roi_mask;
                end
                if nnz(ring) == 0
                    C_neuropil(r,:) = 0;
                else
                    vecn = reshape(ring, [], 1);
                    C_neuropil(r,:) = (vecn' * Yr) / nnz(vecn);
                end
            end
            % subtract scaled neuropil
            C_raw_ns = C_raw - NEUROPIL_SCALE * C_neuropil;
            % Prevent negative fluorescence after subtraction (clip to small positive)
            C_raw_ns(C_raw_ns < 0) = 0;
            C_raw_used = C_raw_ns;
        else
            C_neuropil = [];
            C_raw_used = C_raw;
        end
        
        % ---------- Compute dF/F ----------
        % Df: baseline value(s). If DF_WINDOW provided, compute running percentile baseline per neuron.
        if isempty(DF_WINDOW) || DF_WINDOW <= 0 || DF_WINDOW >= Tframes
            % global percentile per ROI
            Df = prctile(C_raw_used, DF_PRCTILE, 2);    % K x 1
            % avoid zeros
            Df(Df <= 0) = eps;
            C_df = bsxfun(@rdivide, bsxfun(@minus, C_raw_used, Df), Df);  % (K x T)
        else
            % running percentile baseline: compute per-ROI using sliding windows
            Df_mat = zeros(numROIs, Tframes);
            halfwin = floor(DF_WINDOW/2);
            for r = 1:numROIs
                trace = C_raw_used(r,:);
                % pad edges by repeating end values
                trace_p = padarray(trace, [0 halfwin], 'replicate', 'both');
                df_run = zeros(1, Tframes);
                for t = 1:Tframes
                    win = trace_p(t:(t + DF_WINDOW - 1));
                    df_run(t) = prctile(win, DF_PRCTILE);
                end
                Df_mat(r,:) = df_run;
            end
            % avoid zeros
            Df_mat(Df_mat <= 0) = eps;
            C_df = (C_raw_used - Df_mat) ./ Df_mat;
            Df = Df_mat; % store full matrix
        end
        
        % ---------- Prepare CNMF-style outputs for compatibility ----------
        % A: d x K sparse spatial masks (binary)
        A = A_sparse;
        % C_raw: K x T (as computed)
        C_raw_out = C_raw_used;
        % Cn: summary image (use filtered max projection normalized)
        Cn = mat2gray(max_red);
        % Df (baseline) and C_df already computed
        
        % --------------------------------------------
        % Deconvolution (CNMF-E style)
        % --------------------------------------------
        % Use same options structure you already have for green channel
        deconv_options = neuron.options.deconv_options;  
        [C_dec, S_dec, C_denoised, kernel_pars, noise_sn] = ...
            deconv_manual(C_raw_out, deconv_options, true, 'psd');

        % Replace C with deconvolved version
        C = C_dec;
        S = S_dec;
        C_raw_dn = C_denoised;
        
        % ---------- Save results (two files) ----------
        if ~exist(proc_dir,'dir'), mkdir(proc_dir); end
        
        % 1) user-friendly manual results (keeps polygon coords and masks)
        Red_ROIs_save = Red_ROIs;
        Red_Masks_save = Red_Masks;
        red_traces = C_raw_out;
        dFF_red = C_df;
        savefull = fullfile(proc_dir, 'RED_manual_results.mat');
        save(savefull, 'Red_ROIs_save', 'Red_Masks_save', 'red_traces', 'dFF_red', 'Red_MC', '-v7.3');
        fprintf('Saved RED manual results: %s\n', savefull);
        
        % 2) CNMF-style file (for downstream compatibility)
        CNMFe_style_path = fullfile(proc_dir, 'RED_manual_CNMFE_style.mat');
        % save A as sparse, C_raw, C, S, Cn, Df, C_df
        save(CNMFe_style_path, ...
            'A', 'C', 'C_raw_out', 'C_raw_dn', 'S', ...
            'Cn', 'Df', 'C_df', 'kernel_pars', 'noise_sn', ...
            '-v7.3');
        fprintf('Saved RED CNMF-style results: %s\n', CNMFe_style_path);
        
        % ---------- Save outlines on top of max projection ----------
        if SAVE_OUTLINE_PNG
            hf = figure('Name', 'Outlines RED', 'NumberTitle','off', 'Visible','off');
            imagesc(Cn); axis image off; colormap bone; hold on;
            for r = 1:numROIs
                mask = Red_Masks{r};
                B = bwboundaries(mask);
                for b = 1:numel(B)
                    plot(B{b}(:,2), B{b}(:,1), 'r', 'LineWidth', 1.2);
                end
            end
            title(['Outlines RED: ' red_name]);
            drawnow;
            saveas(hf, fullfile(proc_dir, 'Outlines_Neurons_RED.png'));
            close(hf);
            fprintf('Saved outline PNG: %s\n', fullfile(proc_dir, 'Outlines_Neurons_RED.png'));
        end
        
        % ---------- cleanup ----------
        clear Red_MC Yr A_sparse Red_ROIs Red_Masks red_traces dFF_red C_raw_out C_df Df C_neuropil max_red
        
    end
end

disp('Finished extracting red ROIs/traces.');
