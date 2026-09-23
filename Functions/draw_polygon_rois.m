function ROIs = draw_polygon_rois(img, clim)
%DRAW_POLYGON_ROIS  Interactively draw an arbitrary number of polygon ROIs
%   (each with an arbitrary number of vertices) on a background image.
%
%   ROIs = draw_polygon_rois(img, clim)
%
%   img   - [H x W] background image to display (e.g. a max projection)
%   clim  - optional [lo hi] display range, e.g. from prctile(img(:),[1 99])
%
%   Controls: click to place each vertex of the current polygon,
%   double-click (or press Enter) to close it. You'll then be asked
%   whether to draw another. You can also delete the most recently drawn
%   ROI if you made a mistake.
%
%   Returns a struct array with fields:
%     ID       - integer neuron ID (1, 2, 3, ...)
%     Vertices - [n x 2] polygon vertices as drawn (x, y)
%     Mask     - [H x W] logical mask
%     Centroid - [x y] centroid of the mask

    H = size(img, 1);
    W = size(img, 2);

    fig = figure('Name', 'Draw ROIs', 'NumberTitle', 'off');
    if nargin >= 2 && ~isempty(clim)
        imagesc(img, clim);
    else
        imagesc(img);
    end
    axis image; colormap gray; hold on;
    title('Click to place vertices, double-click to close the polygon.');

    ROIs = struct('ID', {}, 'Vertices', {}, 'Mask', {}, 'Centroid', {});
    roi_count = 0;

    while true
        try
            h = drawpolygon('Color', 'y', 'LineWidth', 1.5);
        catch
            break   % figure closed or drawing cancelled
        end

        if ~isvalid(h) || isempty(h.Position) || size(h.Position, 1) < 3
            break   % user cancelled without completing a polygon
        end

        roi_count = roi_count + 1;
        mask = poly2mask(h.Position(:,1), h.Position(:,2), H, W);
        stats = regionprops(mask, 'Centroid');

        ROIs(roi_count).ID = roi_count;
        ROIs(roi_count).Vertices = h.Position;
        ROIs(roi_count).Mask = mask;
        ROIs(roi_count).Centroid = stats(1).Centroid;

        text(stats(1).Centroid(1), stats(1).Centroid(2), num2str(roi_count), ...
             'Color', 'g', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

        choice = questdlg( ...
            sprintf('ROI %d saved. Draw another?', roi_count), ...
            'Continue?', 'Draw another', 'Undo last ROI', 'Finished', 'Draw another');

        switch choice
            case 'Undo last ROI'
                ROIs(roi_count) = [];
                roi_count = roi_count - 1;
                % redraw the figure to remove the stray label/outline
                cla;
                if nargin >= 2 && ~isempty(clim)
                    imagesc(img, clim);
                else
                    imagesc(img);
                end
                axis image; colormap gray; hold on;
                for k = 1:roi_count
                    plot([ROIs(k).Vertices(:,1); ROIs(k).Vertices(1,1)], ...
                         [ROIs(k).Vertices(:,2); ROIs(k).Vertices(1,2)], 'y-', 'LineWidth', 1.5);
                    text(ROIs(k).Centroid(1), ROIs(k).Centroid(2), num2str(k), ...
                         'Color', 'g', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                end
            case 'Finished'
                close(fig);
                return
            otherwise
                % 'Draw another' or dialog closed: loop continues
        end

        if ~isvalid(fig)
            break
        end
    end

    if isvalid(fig)
        close(fig);
    end

    fprintf('Finished: %d ROIs drawn.\n', numel(ROIs));
end

