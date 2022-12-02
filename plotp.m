function plotp(mesh, names, pvec, cmap)

% Ugly function hard coded for M2 and M4 (IAHR Abstract)
figure('units','normalized','outerposition',[0 0 1 1])
np = length(pvec)/mesh.ncells;
p = reshape(pvec, [np, mesh.ncells])';
t = tiledlayout(6,10);
% sgtitle([U.BN, ' : ', ti])
for i = 1:prod(t.GridSize)
    nexttile(i);
    var = p(:,i);
    mesh.plot(var)
    title(names{i})
    colormap(gca, cmap)
    amax = max(abs(var(:,1)), [], 'omitnan') + 1e-3;
    caxis([-amax, amax])
    colorbar;
    axis tight
    %     title(sprintf('%s, %s, %s', names{i}, names{i+np(1)}, names{i + np(1) + np(2)}))
    set(gca, 'XDir','reverse') % Very important
end
end

