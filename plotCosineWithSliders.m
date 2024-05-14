function plotCosineWithSliders()

% Initial values for amplitude and phase
amplitude = 1;
phase = 0;

% Create a figure
fig = figure('Name', 'Cosine Function with Sliders', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);

% Create sliders for amplitude and phase
ampSlider = uicontrol('Style', 'slider', 'Min', 0, 'Max', 5, 'Value', amplitude, 'Position', [20, 20, 200, 20], 'Callback', @updatePlot);
phaseSlider = uicontrol('Style', 'slider', 'Min', -pi, 'Max', pi, 'Value', phase, 'Position', [240, 20, 200, 20], 'Callback', @updatePlot);

% Create axes for plotting
ax = axes('Parent', fig, 'Position', [0.1, 0.3, 0.8, 0.6]);

% Plot the initial cosine function
x = linspace(0, 4 * pi, 1000);
y = amplitude * cos(x - phase);
plot(ax, x, y, 'LineWidth', 2);
title('Cosine Function with Sliders');
xlabel('Time');
ylabel('Amplitude');

    function updatePlot(~, ~)
        % Get current slider values
        amplitude = get(ampSlider, 'Value');
        phase = get(phaseSlider, 'Value');

        % Update the plot
        y = amplitude * cos(x - phase);
        plot(ax, x, y, 'LineWidth', 2);

        % Update title
        title(ax, sprintf('Cosine Function with Sliders\nAmplitude: %.2f, Phase: %.2f', amplitude, phase));
    end



end
