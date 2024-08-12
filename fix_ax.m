function ax = fix_ax(ax, opts)
%ax = nexttile;
set(ax, 'FontName',  "BookAntiqua")
set(ax, 'FontSize',  10)
axis tight
ax.TickLabelInterpreter = 'latex';
if strcmp(opts.ylim, "sym")
    YL = get(gca, 'YLim');
    maxlim = max(abs(YL));
    set(gca, 'YLim', [-maxlim maxlim]);
elseif strcmp(opts.ylim, "pos")
    YL = get(gca, 'YLim');
    maxlim = max(abs(YL));
    set(gca, 'YLim', [0 maxlim]);
end
grid on
end