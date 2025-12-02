%% Manually extract cells from RED movie
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
        Red_MC = fullfile(proc_dir, 'MC_red.mat');

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

            for tt = 1:Tframes
                frame = Red_MC(:,:,tt);           % extract frame first
                red_traces(r,tt) = mean(double(frame(idxPix)));
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

disp('Finished extracting red ROIs/traces.');