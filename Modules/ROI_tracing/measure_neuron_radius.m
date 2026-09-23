function radius_px = measure_neuron_radius(img, clim)
%MEASURE_NEURON_RADIUS  Interactively measure a typical neuron radius (in
%   pixels) by clicking across the diameter of a few clearly visible
%   somas. Use this on the channel/image where cells are easiest to see
%   unambiguously (usually the green channel, given well-visible
%   pyramidal neurons) - the resulting pixel radius reflects the optics
%   and resolution of the setup, so the same value is a reasonable
%   starting point for tuning Bg_Sigma on the red channel too.
%
%   radius_px = measure_neuron_radius(img, clim)
%
%   img  - [H x W] background image to click on, e.g. mean(MC_f.MC(:,:,1:200),3)
%   clim - optional [lo hi] display range, e.g. from prctile(img(:),[1 99])
%
%   Controls: for each neuron, click twice - once on each side of the
%   soma, spanning its diameter. Repeat for several (5-10 recommended)
%   clearly visible, roughly average-sized neurons. Press Enter (without
%   clicking) when done.
%
%   Returns radius_px, the average of all measured radii.

    figure('Name', 'Measure neuron radius');
    if nargin >= 2 && ~isempty(clim)
        imagesc(img, clim);
    else
        imagesc(img);
    end
    axis image; colormap gray; hold on;
    title('Click 2 points across each neuron''s diameter. Press Enter when done.');

    diameters = [];
    while true
        [x, y] = ginput(2);
        if numel(x) < 2
            break   % user pressed Enter without completing a pair
        end
        plot(x, y, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
        d = sqrt(diff(x)^2 + diff(y)^2);
        diameters(end+1) = d; %#ok<AGROW>
        fprintf('Measurement %d: diameter = %.1f px (radius %.1f px)\n', numel(diameters), d, d/2);
    end

    if isempty(diameters)
        warning('No measurements taken. Returning NaN.');
        radius_px = NaN;
        return
    end

    radius_px = mean(diameters) / 2;
    fprintf('\nAverage over %d measurements: radius = %.1f px (std %.1f px)\n', ...
            numel(diameters), radius_px, std(diameters)/2);
end