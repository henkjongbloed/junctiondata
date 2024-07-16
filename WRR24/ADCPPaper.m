%% ADCP data: Method paper Jongbloed et al. 2024

% 1. Solution needs Solver as input
% 2. Fix plotting function, with possibility of amp-phase
% 3. Implement K-fold CV
% 4. Implement sens. to noise
% 5. K-fold CV for 3 different meshes, 3 diff points in 2D par. space
% 6. Hydrostatic press using decomposition function, not on mesh.
% 7. Rewrite statistical section

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextFontSize', 10);
%set(groot, 'default

clearvars % Clear all variables from the workspace
close all
figsave = 0;

RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'WP2\TwoJunctions\code'])

tak = 1;

DS = 'OMHA14';

%Import preprocessed semi-raw data (ADCP / CTD / H structs)
[adcp, ctd, h] = import_data(RF, DS);
V = rdi.VMADCP(adcp{1});

%% Geometry + Time
constituents = {'M2', 'M4'};
%
water_level = VaryingWaterLevel(datetime(datevec(h.t)), h.wl); % set water level to varying
water_level.model = TidalScalarModel(constituents = constituents);
water_level.model.scalar_name = 'eta'; % Scalar
water_level.get_parameters();
V.water_level_object = water_level;
load(strcat("xs_", DS)); % define the cross-section
load(strcat("ef_", DS));
xs = xs(tak);
ef = ef(tak); % -> define further fitering such that nensembles equals number of ef ensembles


V.filters = Filter;
bathy = BathymetryScatteredPoints(ef, V);
interpolators = [bathy.interpolator];
[interpolators.span] = deal(0.01);

mesh_maker = SigmaZetaMeshFromVMADCP(ef, xs, bathy, 'NoExpand', V);

%mesh_mean = mesh_maker.get_mesh(resn = 50, resz = 20); % get mesh at mean water level
mesh_mean = mesh_maker.get_mesh(resn = 25 , resz = 13);
%% Plot 1: Mesh

%% SolverOptions
opts = SolverOptions(extrapolate_vert = 0);

%% Empirical model: TaylorTidal
flow_model = TaylorTidalVelocityModel; % VelocityModel, thus 3 Cartesian components
flow_model.constituents = constituents;
flow_model.s_order = [1 1 1];
flow_model.n_order = [1 1 1];
flow_model.sigma_order = [1 1 1];

%% Solver options and regularization
flow_regs = regularization.Velocity.get_all_regs(mesh_mean, bathy, xs, flow_model, opts, 'NoExpand', V);

flow_regs(1).weight = [0, 30, 250];
flow_regs(2).weight = [0, 30, 250];
flow_regs(3).weight = [0, 3, 250];
flow_regs(4).weight = [0, 3, 250];
flow_regs(5).weight = [0, 30, 250];

%% Solve and plot

flow_solv = LocationBasedVelocitySolver(mesh_mean, bathy, xs, ef, flow_model, opts, 'NoExpand', V, flow_regs);
flow_solv.rotation = xs.angle;
flow = flow_solv.get_solution();

%% Plot state
names = {flow_model.all_names{[1,2,3,4,5,16,21,41,6,31,56]}};

% names = {flow_model.all_names{[1,2,3,22,42,7,32,57]}};

flow.plot_solution({flow_model.all_names{[1,2,3,4,5,21,41]}}, representation = "Aphi", sol_idx = 1:3) % u,v,w

flow.plot_solution({flow_model.all_names{15+[1,2,4]}}, representation = "Aphi", sol_idx = 1:3)

flow.plot_solution({flow_model.all_names{[16,17,19, 6,31,56]}}, representation = "Aphi", sol_idx = 1:3)
%% Plot CV as function of l_s
nreg = 30;
lc = 5;
ls = logspace(-4,3,nreg-1)';
regp = [[0;lc*ones(nreg-1,1)], 0*ones(nreg,1), [0;ls], [0;ls], [0;lc*ones(nreg-1,1)]];
clear regpc
for i =1:nreg
    regpc{i} = regp(i,:);
end

CV = flow.cross_validate_single(regpc);

save("cvfine", 'CV')

semilogx([0;ls], [CV{1:end,1}]/CV{1,1})
hold on
semilogx(ls, [CV{2:end,2}]/CV{1,2})



% for i = 1:size(flow.p,2)
%     figure;
%     mesh_mean.plot(flow.pars(:,1,i)); colorbar;
%     mesh_mean.plot(flow.pars(:,2,i)); colorbar;
%     mesh_mean.plot(flow.pars(:,3,i)); colorbar;
%     mesh_mean.plot(flow.pars(:,4,i)); colorbar;
%     mesh_mean.plot(flow.pars(:,5,i)); colorbar;
%
%     figure;
%     mesh_mean.plot(flow.pars(:,21,i)); colorbar;
%     mesh_mean.plot(flow.pars(:,22,i)); colorbar;
%     mesh_mean.plot(flow.pars(:,23,i)); colorbar;
%     mesh_mean.plot(flow.pars(:,24,i)); colorbar;
%     mesh_mean.plot(flow.pars(:,25,i)); colorbar;
%
%     figure;
%     mesh_mean.plot(flow.pars(:,41,i)); colorbar;
%     mesh_mean.plot(flow.pars(:,42,i)); colorbar;
%     mesh_mean.plot(flow.pars(:,43,i)); colorbar;
%     mesh_mean.plot(flow.pars(:,44,i)); colorbar;
%     mesh_mean.plot(flow.pars(:,45,i)); colorbar;
% end
%% Mesh Plot
m = makefigure(20,3);
t = tiledlayout('flow', TileSpacing="tight", Padding="tight");
ax=nexttile;
%subplot(1,2,1)
flow.plot_mesh(ax, vertical = 'z');
title('$z$ - coordinates', 'Interpreter','latex','FontSize', 10);
xlabel("y")
ax=nexttile;
%subplot(1,2,2)
flow.plot_mesh(ax, vertical = 'sig');
title('$\sigma$ - coordinates', 'Interpreter','latex','FontSize', 10);
xlabel("y")


%ax=nexttile;

%gca.TickLabelInterpreter = 'latex';
%% Mesh Plot
m = makefigure(2*16/3,2*3);
t = tiledlayout('flow', TileSpacing="tight", Padding="tight");
% ax=nexttile;
% %subplot(1,2,1)
% flow.plot_mesh(ax, vertical = 'z');
% title('$z$ - coordinates', 'Interpreter','latex','FontSize', 10);
% %xlabel([])
% ax=nexttile;
% %subplot(1,2,2)
% flow.plot_mesh(ax, vertical = 'sig');
% title('$\sigma$ - coordinates: ', 'Interpreter','latex','FontSize', 10);
% %xlabel([])
ax=nexttile;
%subplot(1,2,2)
flow.plot_mesh(ax, vertical = 'sig');
title('$\lambda_{c}$, $\lambda_{s}$ , $\lambda_{N}$ , $\lambda_{L}$ , $\lambda_{H}$', 'Interpreter','latex','FontSize', 10);
%xlabel([])
%gca.TickLabelInterpreter = 'latex';
%% Mesh Plot - use to obtain better figs!
m = makefigure(2*16/3,2*3);
t = tiledlayout('flow', TileSpacing="tight", Padding="tight");
% ax=nexttile;
% %subplot(1,2,1)
% flow.plot_mesh(ax, vertical = 'z');
% title('$z$ - coordinates', 'Interpreter','latex','FontSize', 10);
% %xlabel([])
% ax=nexttile;
% %subplot(1,2,2)
% flow.plot_mesh(ax, vertical = 'sig');
% title('$\sigma$ - coordinates: ', 'Interpreter','latex','FontSize', 10);
% %xlabel([])
ax=nexttile;
%subplot(1,2,2)
flow.plot_mesh(ax, vertical = 'sig');
title('$ = $     $E^{train}_b(\mathbf{\lambda}) / E^{train}_b(\mathbf{0})$ , $E^{gen}_b(\mathbf{\lambda}) / E^{gen}_b(\mathbf{0})$ , $\log_{10}E_p(\mathbf{\lambda})$', 'Interpreter','latex','FontSize', 10);
%xlabel([])
%gca.TickLabelInterpreter = 'latex';

%% Mesh Plot - use to obtain better figs!
m = makefigure(2*16/3,2*3);
t = tiledlayout('flow', TileSpacing="tight", Padding="tight");
% ax=nexttile;
% %subplot(1,2,1)
% flow.plot_mesh(ax, vertical = 'z');
% title('$z$ - coordinates', 'Interpreter','latex','FontSize', 10);
% %xlabel([])
% ax=nexttile;
% %subplot(1,2,2)
% flow.plot_mesh(ax, vertical = 'sig');
% title('$\sigma$ - coordinates: ', 'Interpreter','latex','FontSize', 10);
% %xlabel([])
ax=nexttile;
%subplot(1,2,2)
flow.plot_mesh(ax, vertical = 'sig');
title('MSE     $E^{train}_b(\mathbf{\lambda}) / E^{train}_b(\mathbf{0})$ , $E^{gen}_b(\mathbf{\lambda}) / E^{gen}_b(\mathbf{0})$ , $\log_{10}E_p(\mathbf{\lambda})$', 'Interpreter','latex','FontSize', 10);
%xlabel([])
%gca.TickLabelInterpreter = 'latex';

%% Evaluate to visualize Dw/Dt and hydrostatic press.
m = makefigure(20, 4);
t = tiledlayout(m, 1, 3, TileSpacing = "tight", TileIndexing = "columnmajor");
lam = {'\mathbf{\lambda}_N', '\mathbf{\lambda}_L', '\mathbf{\lambda}_H'};
% titles =
for row = 1:3
    ax=nexttile;
    [acc,X] = get_acc(flow,row, .5);
    Dsum = sum(cat(4, acc{:}), 4);
    ncolor = 100;
    amax = max(abs([min(Dsum, [], "all"), max(Dsum, [], "all")]));
    levels = linspace(-amax, amax, ncolor);
    [~,ha]=contourf(squeeze(X.Y(1,:,:))', squeeze(X.Z(1,:,:))', Dsum' , "LineColor",'none');
    c = colorbar;
    c.TickLabelInterpreter = 'latex';
    c.FontSize = 10;
    title(['$Dw/Dt: \hspace{.1cm} \lambda = ', lam{row}, '$'], 'interpreter', 'latex', 'FontSize', 12);
    colormap(ax, helpers.cmaps("velmap"))
    clim([-amax, amax])
    ylim([min(X.Z, [], 'all'), max(X.Z, [], 'all')])
    set(gca, 'XDir','reverse') % Very important

    set(gca, 'XTick',[])
    set(gca, 'YTick',[])
    %set(ax, "YLabel", "z")
    %set(ax, "XLabel", "y")
end
t.XLabel.String = 'y';
t.YLabel.String = 'z';
t.XLabel.Interpreter = 'latex';
t.YLabel.Interpreter = 'latex';
