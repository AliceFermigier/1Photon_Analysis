%% ManualROIsRED.m
% Manual ROI definition for RED channel using repository functions:
%  - ROI_hand_drawn.m
%  - extract_traces.m
%  - compute_dff.m  (if available; otherwise we implement the repo version)

if ~exist('All_nam','var')
    error('All_nam must be in workspace (same one used during motion correction).');
end

for p = 1:size(All_nam,1)
    subfolders = All_nam{p};

    for pp = 1:numel(subfolders)
        Current_folder = [subfolders{pp} filesep];
        proc_dir = fullfile(Current_folder, 'processed_data');
        red_mat_path = fullfile(proc_dir, 'MC_red.mat');

        if ~exist(red_mat_path,'file')
            warning('MC_red.mat not found in %s. Skipping.\n', Current_folder);
            continue
        end

        % Load corrected movie
        S = load(red_mat_path);
        if ~isfield(S,'Red_MC')
            error('MC_red.mat found but does not contain Red_MC.');
        end
        Y = double(S.Red_MC);          % d1 x d2 x T
        [d1, d2, T] = size(Y);
        d = d1*d2;

        % ----- Create summary image for drawing -----
        % Use same approach as green channel workflow
        summary_img = max(Y, [], 3);

        figure; imagesc(summary_img); axis image off; colormap grey;
        title('Draw RED-Channel ROIs (double click to finish ROI)');
        drawnow;

        % ----- Use ROI_hand_drawn from repo -----
        % This function returns masks as a sparse A matrix (d x K)
        fprintf('Draw ROIs, close figure when finished.\n');
        A_red = ROI_hand_drawn(summary_img);   % repo function
        K = size(A_red,2);
        close(gcf);

        % ----- Extract fluorescence traces using repo function -----
        % extract_traces expects Y reshaped as d x T
        Yr = reshape(Y, d, T);

        % returns raw fluorescence C_raw_red (K x T)
        C_raw_red = extract_traces(Yr, A_red);

        % ----- Compute ?F/F using repository standard algorithm -----
        % fallback to repository-consistent dF/F
        F0_red = prctile(C_raw_red, 20, 2);
        C_dff_red = (C_raw_red - F0_red) ./ F0_red;

        % Save everything in repo-standard variables
        save(fullfile(proc_dir,'RED_manual_results.mat'), ...
             'A_red','C_raw_red','C_dff_red','F0_red','summary_img','-v7.3');

        fprintf('Saved RED manual ROI results ? %s\n', ...
                 fullfile(proc_dir,'RED_manual_results.mat'));

        clear Y Yr A_red C_raw_red C_dff_red F0_red
    end
end

disp('Finished manual RED ROI extraction.');
