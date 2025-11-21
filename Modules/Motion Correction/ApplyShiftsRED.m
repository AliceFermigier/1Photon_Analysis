%% Apply saved green motion-correction shifts to RED movie + manual ROI extraction
% Usage: run after Motion_correction has finished for the session(s).
% Expects each session folder to contain:
%   processed_data/MC_Shifts.mat   (saved by Motion_correction)
% and a raw red movie file named like: 839R_EPM.tif or 839R_EPM.tiff
%
% Globals expected from your pipeline:
%   All_nam            cell array from uipickfiles, same shape you used before
%   Data_Format_in     'tif' or 'tiff' etc.
%   T_DS_factor        temporal downsampling factor (same as used for green MC)
%   Spatial_Downsampling (or Resize_factor)   spatial resize factor used for MC
%
% Example: set in base workspace then run this script.
%   Spatial_Downsampling = 0.25;
%   T_DS_factor = 1;
%   ApplyShiftsAndRedROIs

if ~exist('All_nam','var')
    error('All_nam is not in workspace. Provide same All_nam used for motion correction.');
end
if ~exist('Data_Format_in','var') || ~exist('T_DS_factor','var') || ~exist('Spatial_Downsampling','var')
    error('Make sure Data_Format_in, T_DS_factor and Spatial_Downsampling are defined in the workspace.');
end

Resize_factor = Spatial_Downsampling; % match Motion_correction arg name

% Loop over animals / sessions
for p = 1:size(All_nam,1)
    subfolders = All_nam{p}; % cell array of session folder paths for this animal
    
    for pp = 1:numel(subfolders)
        Current_folder = [subfolders{pp} filesep];
        proc_dir = fullfile(Current_folder, 'processed_data');
        shifts_file = fullfile(proc_dir, 'MC_Shifts.mat');
        if ~exist(shifts_file,'file')
            warning('No MC_Shifts.mat in %s — skipping this session.', Current_folder);
            continue
        end
        
        % Load the shift collection saved by Motion_correction
        S = load(shifts_file, 'Shift_collection');
        Shift_collection = S.Shift_collection; clear S;
        
        % Find the RED raw movie file
        listing = dir(fullfile(Current_folder, ['*.' Data_Format_in]));
        names = {listing.name};
        % Match pattern like 839R_EPM.tif  (capital or small R)
        red_mask = ~cellfun('isempty', regexp(names, '^\d+R_.*\.tif{1,2}$', 'once'));
        if ~any(red_mask)
            % fallback: any file containing 'R' by convention
            red_mask = contains(names, 'R_') | contains(names, 'r_');
        end
        if ~any(red_mask)
            warning('No red file found in %s. Files: %s. Skipping.', Current_folder, strjoin(names,','));
            continue
        end
        red_name = names{find(red_mask,1,'first')};
        red_path = fullfile(Current_folder, red_name);
        fprintf('Processing RED movie: %s\n', red_path);
        
        % --------- Load RED movie and resize (same as in Motion_correction) ----------
        switch lower(Data_Format_in)
            case {'tif','tiff'}
                % Use loadtiff (your pipeline has it). This returns HxWxT
                raw_red = loadtiff(red_path);   % may be uint16 or uint8
                % resize spatially exactly like Motion_correction
                sz = size(raw_red);
                try
                    % If original code used imresize on whole stack, do same:
                    red_resized = imresize(raw_red, [sz(1)*Resize_factor, sz(2)*Resize_factor], 'bicubic');
                catch
                    % If memory heavy, resize frame-by-frame
                    Tframes = size(raw_red,3);
                    red_resized = zeros(round(sz(1)*Resize_factor), round(sz(2)*Resize_factor), Tframes, class(raw_red));
                    for tt = 1:Tframes
                        red_resized(:,:,tt) = imresize(raw_red(:,:,tt), [round(sz(1)*Resize_factor), round(sz(2)*Resize_factor)], 'bicubic');
                    end
                end
                
            case 'hdf5'
                % If you used h5 in your pipeline you will need to adapt this block
                data = h5read(red_path, '/images');
                red_resized = imresize(data, [size(data,1)*Resize_factor, size(data,2)*Resize_factor], 'bicubic');
                clear data;
                
            case 'mat'
                tmp = matfile(red_path);
                data = tmp.data;
                red_resized = imresize(data, [size(data,1)*Resize_factor, size(data,2)*Resize_factor], 'bicubic');
                clear data;
                
            otherwise
                error('Data format %s not implemented in this helper script.', Data_Format_in);
        end
        
        % Convert to double in range 0..1 in same fashion as Motion_correction did:
        % Motion_correction used mat2gray(A, [0 Max_Value]) where Max_Value depended on class
        switch class(red_resized)
            case 'uint16'
                Max_Value = 2^16 - 1;
            case 'uint8'
                Max_Value = 2^8 - 1;
            otherwise
                Max_Value = double(max(red_resized(:)));
        end
        red_resized = mat2gray(red_resized, [0 Max_Value]);  % now double 0..1
        
        % Temporal downsampling: keep the same frames as MRDFT (1:T_DS_factor:end)
        if T_DS_factor > 1
            red_ds = red_resized(:,:,1:T_DS_factor:end);
        else
            red_ds = red_resized;
        end
        clear red_resized raw_red;
        
        % Preallocate corrected movie array (same class used in MC_f.MC -> uint16)
        Red_MC = zeros(size(red_ds), 'like', im2uint16(red_ds));
        
        % Determine chunk boundaries used by MC to iterate correctly
        % Motion_correction created Shift_collection with length = numFiles (chunks)
        numChunks = numel(Shift_collection);
        
        % We must know how many frames are in each chunk in the downsampled movie.
        % The Shift_collection{c} is an N x 2 matrix where N equals number of frames in that chunk after downsampling.
        frameStart = 1;
        for c = 1:numChunks
            Shifts_all = Shift_collection{c}; % nFrames_chunk x 2
            nF = size(Shifts_all,1);
            if nF == 0
                continue
            end
            idxFrames = frameStart:(frameStart + nF - 1);
            frameStart = frameStart + nF;
            % safety check
            if idxFrames(end) > size(red_ds,3)
                warning('Chunk frames exceed the downsampled red movie length in %s. Truncating.', Current_folder);
                idxFrames = idxFrames(idxFrames <= size(red_ds,3));
                nF = numel(idxFrames);
                Shifts_all = Shifts_all(1:nF,:);
            end
            
            % Apply per-frame translations using same shifts
            for ff = 1:nF
                frameIdx = idxFrames(ff);
                thisShift = Shifts_all(ff,:); % [y x] as in Image_Registration_gray
                % imtranslate expects [x y] translation vector when using 'XData','YData'?
                % Using imtranslate(I, [tx ty]) where vector is [x y]
                % Our shift is [row_shift, col_shift] -> [y x], so swap order
                tx = thisShift(2);
                ty = thisShift(1);
                % apply translation; use 'FillValues' 0 to avoid NaNs
                corrected = imtranslate(red_ds(:,:,frameIdx), [tx, ty], 'bicubic', 'FillValues', 0);
                Red_MC(:,:,frameIdx) = im2uint16(corrected);
            end
        end
        
        % Save corrected red movie to processed_data
        if ~exist(proc_dir,'dir'), mkdir(proc_dir); end
        savefast(fullfile(proc_dir, 'MC_red.mat'), 'Red_MC');  % large file
        fprintf('Saved corrected red movie: %s\n', fullfile(proc_dir, 'MC_red.mat'));
        
        % Generate mean projection for ROI drawing (use the uint16 -> double)
        mean_red = mean(double(Red_MC), 3);
        figure; imagesc(mean_red); axis image; colormap bone; title(['Mean RED: ' red_name]);
        
        % Manual ROI drawing (polygons). Store ROI masks and positions
        Red_ROIs = {};
        Red_Masks = {};
        keep_drawing = true;
        uiwait(msgbox('Draw ROIs on the RED mean projection. Double-click to finish each ROI. Click OK to start.', 'Draw ROIs','modal'));
        while keep_drawing
            h = drawpolygon('LineWidth',1.5,'Color','r');
            pos = h.Position;
            Red_ROIs{end+1} = pos;
            mask = poly2mask(pos(:,1), pos(:,2), size(Red_MC,1), size(Red_MC,2));
            Red_Masks{end+1} = mask;
            choice = questdlg('Add another ROI?', 'ROI','Yes','No','Yes');
            if strcmp(choice,'No'), keep_drawing = false; end
        end
        close(gcf);
        
        % Extract traces for each ROI from the corrected red movie
        numROIs = numel(Red_Masks);
        Tframes = size(Red_MC,3);
        red_traces = zeros(numROIs, Tframes);
        for r = 1:numROIs
            m = Red_Masks{r};
            idxPix = find(m);
            % extract mean across pixels for each timepoint
            % Red_MC stored as uint16 -> convert to double and normalize to 0..1 using Max_Value
            for tt = 1:Tframes
                red_traces(r,tt) = mean(double(Red_MC(:,:,tt)(idxPix)));
            end
        end
        
        % Normalize traces to 0..1 by dividing by max uint16
        red_traces = red_traces ./ double(intmax('uint16'));
        
        % Compute simple dF/F using 10th percentile as baseline
        F0 = prctile(red_traces, 10, 2);
        dFF_red = (red_traces - F0) ./ F0;
        
        % Save ROIs & traces
        savefull = fullfile(proc_dir, 'RED_manual_results.mat');
        save(savefull, 'Red_ROIs', 'Red_Masks', 'red_traces', 'dFF_red', 'Red_MC', '-v7.3');
        fprintf('Saved RED ROI results to %s\n', savefull);
        
        % Clear large variables before next session
        clear Red_MC red_ds red_traces dFF_red Red_ROIs Red_Masks mean_red
    end
end

disp('Finished applying green shifts to red channel and extracting red ROIs/traces.');
