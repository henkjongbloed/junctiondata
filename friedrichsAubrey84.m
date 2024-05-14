%FriedrichsAubrey1984
close all
clearvars

A2 = 1; P2 = 0;

A4 = 0; P4 = 0;

t = 0:(2/60):(2*24);

T = 12.42;
om = 2*pi/T;

y2 = simpleCurve(t, A2, om, P2);
y4 = simpleCurve(t, A4, 2*om, P4);

fig = figure('Name', 'Friedrichs and Aubrey 1988', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
% Create sliders for amplitude and phase
combinedPlot(t, A2, A4, om, P2, P4)
PSlider = uicontrol('Style', 'slider', 'Min', -pi, 'Max', pi, 'Value', P4, 'Position', [240, 20, 200, 20], 'Callback', @(src, event) updateP(src, event, ASlider,  t, A2, om, P2));
ASlider = uicontrol('Style', 'slider', 'Min', 0, 'Max', 2, 'Value', A4, 'Position', [20, 20, 200, 20], 'Callback', @(src, event) updateA(src, event, PSlider, t, A2, om, P2));
PSlider = uicontrol('Style', 'slider', 'Min', -pi, 'Max', pi, 'Value', P4, 'Position', [240, 20, 200, 20], 'Callback', @(src, event) updateP(src, event, ASlider,  t, A2, om, P2));

function updateP(src, event, ASlider, t, A2, om, P2)
    % Get current slider values
    A4 = get(ASlider, 'Value');
    P4 = get(src, 'Value');

        % Update the plot
    combinedPlot(t, A2, A4, om, P2, P4)


    % Update title
%     title(ax, sprintf('Cosine Function with Sliders\nAmplitude: %.2f, Phase: %.2f', amplitude, phase));
end

function updateA(src, event, PSlider, t, A2, om, P2)
    % Get current slider values
    A4 = get(src, 'Value');
    P4 = get(PSlider, 'Value');

        % Update the plot
    combinedPlot(t, A2, A4, om, P2, P4)


    % Update title
%     title(ax, sprintf('Cosine Function with Sliders\nAmplitude: %.2f, Phase: %.2f', amplitude, phase));
end



function dat = simpleCurve(t, A, om, P)
    dat = A.*cos(om.*t - P);
    %plot(ax, t, dat, 'LineWidth', 2);
end

function combinedPlot(t, A2, A4, om, P2, P4)
    ax = gca;
    for i = 1:numel(ax.Children)
        ax.Children(i).XData = [];
        ax.Children(i).YData = [];
    end
    hold on
    y2 = simpleCurve(t, A2, om, P2);
    y4 = simpleCurve(t, A4, 2*om, P4);
    plot(t, y2, 'k:')
    plot(t, y4, 'k--')
    plot(t, y2 + y4,'b-', 'LineWidth', 2)
    hold off
    ax.YLim = [-2,2];
    title(['FA88: ', '2\phi_{M2} - \phi_{M4} = ', num2str(2*P2 - P4), ' and A_{M4}/A_{M2} = ', num2str(A4/A2)])
end
