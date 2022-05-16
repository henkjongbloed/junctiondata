function plotParsx(U)

% Ugly function hard coded for M2 and M4 (IAHR Abstract)
figure;
npars = U.T.velocity_model.npars  ;
% ntak = numel(U);
idx_0 = [1, npars(1) + 1, sum(npars(1:2)) + 1];
idx_ampl = [(idx_0(1)+1):2:(idx_0(2)-2), (idx_0(2)+1):2:(idx_0(3)-2), (idx_0(3)+1):2:sum(npars)];
idx_phase = idx_ampl + 1;

t = tiledlayout(3,2);

nexttile([1,2])
U.mesh_mean.plot_vec(U.tid_pars(:,idx_0))
amax = max(max(abs(U.tid_pars(:,idx_0)), [], 'omitnan'));
colormap(gca, U.velmap)
caxis([-amax, amax])
colorbar;
axis tight
title([U.BN, ': Subtidal Flow [m/s]'])
set(gca, 'XDir','reverse') % Very important

nexttile;
U.mesh_mean.plot_vec(U.tid_pars(:,idx_ampl(:,1)));
amax = max(max(abs(U.tid_pars(:,idx_ampl(:,1))), [], 'omitnan'));
colormap(gca, U.velmap(10:end,:))
caxis([0, amax])
colorbar;
title('M2 Amplitude [m/s]')
set(gca,'xticklabel',[])
axis tight
set(gca, 'XDir','reverse') % Very important

nexttile;
U.mesh_mean.plot_vec(U.tid_pars(:,idx_phase(:,1)))
% amax = max(max(abs(U.tid_pars(:,idx_phase)), [], 'omitnan'));
colormap(gca, U.phimap)
caxis([-pi, pi])
colorbar;
title('M2 Phase [rad]')
axis tight
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca, 'XDir','reverse') % Very important

nexttile;
U.mesh_mean.plot_vec(U.tid_pars(:,idx_ampl(:,2)))
amax = max(max(abs(U.tid_pars(:,idx_ampl(:,2))), [], 'omitnan'));
colormap(gca, U.velmap(10:end,:))
caxis([0, amax])
colorbar;
title('M4 Amplitude [m/s]')
set(gca, 'XDir','reverse') % Very important
axis tight


nexttile;
U.mesh_mean.plot_vec(U.tid_pars(:,idx_phase(:,2)))
% amax = max(max(abs(U.tid_pars(:,idx_phase)), [], 'omitnan'));
colormap(gca, U.phimap)
caxis([-pi, pi])
colorbar;
title('M4 Phase [rad]')
axis tight
set(gca,'yticklabel',[])
set(gca, 'XDir','reverse') % Very important

t.XLabel.String = 'y [m]';
t.YLabel.String = 'z [m]';
t.TileSpacing = 'tight';
t.Padding = 'tight';

end