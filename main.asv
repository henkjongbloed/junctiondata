%% ADCP data: Method paper
% Created by Henk Jongbloed

clearvars % Clear all variables from the workspace
close all % Close all figures
figsave = 0;

RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'WP2\TwoJunctions\code'])

ds = 1;
tak = 1;

DS = {'OMHA14','NMOMNW15'}; DS = DS{ds};
BN = {{' Hartel Canal'; ' Old Meuse South'; ' Old Meuse North'},...
    {' New Meuse'; ' Old Meuse'; ' Rotterdam Waterway'}}; BN = BN{ds};

BNs = {{'HK2'; 'OMS2'; 'OMN2'},...
    {'NM2'; 'OM2'; 'RW2'}}; BNs = BNs{ds};

%% Import preprocessed semi-raw data (ADCP / CTD / H structs)

[adcp, ctd, h] = import_data(RF, DS);
for t = 1:numel(tak)
     V(t) = rdi.VMADCP(adcp{tak(t)}); 
end

%% Initial inspection
constituents = {'M2', 'M4'};
% 
water_level = VaryingWaterLevel(datetime(datevec(h.t)), h.wl); % set water level to varying
water_level.model = TidalScalarModel(constituents = constituents);
water_level.model.scalar_name = 'eta'; % Scalar
water_level.get_parameters();
[V.water_level_object] = deal(water_level);

%% create bathymetry
%bathy = BathymetryScatteredPoints(V); % create a bathymetry (includes all data)

dat_selection = "saved";
switch dat_selection
    case "manual"
        [ef, xs] = cross_section_selector(V);
    case "default"
        xs = XSection(V); % define the cross-section
        ef = EnsembleFilter(V);
    case "saved" % Save ensemblefilter and cross-section once
        load(strcat("xs_", DS)); % define the cross-section
        load(strcat("ef_", DS));
        xs = xs(tak);
        ef = ef(tak); % -> define further fitering such that nensembles equals number of ef ensembles
    otherwise 
        xs = XSection(V); % define the cross-section
        ef = EnsembleFilter(V);
end

[V.filters] = deal(Filter);
bathy = BathymetryScatteredPoints(ef, V);
interpolators = [bathy.interpolator];
[interpolators.span] = deal(0.01);

mesh_maker = SigmaZetaMeshFromVMADCP(ef, xs, bathy, 'NoExpand', V);

mesh_mean = mesh_maker.get_mesh(resn = 50, resz = 20); % get mesh at mean water level
mesh_mean.plot()


% %% Mesh plot
% % mesh_mean.plot('FixAspectRatio', false)
% 
% %% Visualize mesh
% % m = makefigure(10,3);
% % mesh_mean(1).plot('sigma', false, 'FixAspectRatio', false);
% % title('$z$ - coordinates', 'Interpreter','latex','FontSize', 12);
% % axis tight
% % ax = gca;
% % set(ax, 'XDir','reverse')
% % ax.TickLabelInterpreter = 'latex';
% % axs=axes(m,'visible','off');
% % %m.Title.Visible='on';
% % axs.XLabel.Visible='on';
% % axs.YLabel.Visible='on';
% % %ylabel(m,'yourYLabel');
% % xlabel(axs,'y [m]', 'Interpreter','latex');
% % ylabel(axs, '$z [m]$', 'Interpreter', 'latex')
% % 
% % 
% % m = makefigure(10,3);
% % mesh_mean(1).plot('sigma', true, 'FixAspectRatio', true);
% % title('$\sigma$ - coordinates', 'Interpreter','latex', 'FontSize', 12);
% % axis tight
% % ax = gca;
% % set(ax, 'XDir','reverse')
% % ax.TickLabelInterpreter = 'latex';
% % axs=axes(m,'visible','off');
% % %m.Title.Visible='on';
% % axs.XLabel.Visible='on';
% % axs.YLabel.Visible='on';
% % %ylabel(m,'yourYLabel');
% % xlabel(axs,'y [m]', 'Interpreter','latex');
% % ylabel(axs,'$\sigma$', 'Interpreter', 'latex')
% 
%% SolverOptions
opts = SolverOptions();


%% Empirical model: TaylorTidal
flow_model = TaylorTidalVelocityModel;
flow_model.constituents = constituents;
flow_model.s_order = deal([1 1 1]);
flow_model.n_order = deal([1 1 1]);
flow_model.sigma_order = deal([1 1 1]);

%% Solver options and regularization

% %% Add regularization terms
flow_regs = regularization.Velocity.get_all_regs(mesh_mean, bathy, xs, flow_model, 'NoExpand', V);

flow_regs(1).weight = [10 1000];
flow_regs(2).weight = [10 1000];
flow_regs(3).weight = [1 100];
flow_regs(4).weight = [1 100];
flow_regs(5).weight = [10 1000];
%% Solve

flow_solv = LocationBasedVelocitySolver(mesh_mean, bathy, xs, ef, flow_model, opts, 'NoExpand', V, flow_regs); 
[flow_solv.rotation] = deal(xs.angle);
flow = flow_solv.get_solution();


%% Empirical model: Tidal
flow_model = TidalVelocityModel;
flow_model.constituents = constituents;
% flow_model.s_order = deal([1 1 1]);
% flow_model.n_order = deal([1 1 1]);
% flow_model.sigma_order = deal([1 1 1]);

%% Solver options and regularization

% %% Add regularization terms
flow_regs = regularization.Velocity.get_all_regs(mesh_mean, bathy, xs, flow_model, 'NoExpand', V);

flow_regs(1).weight = [10 1000];
flow_regs(2).weight = [10 1000];
flow_regs(3).weight = [1 100];
flow_regs(4).weight = [1 100];
flow_regs(5).weight = [10 1000];
%% Solve

flow_solv = LocationBasedVelocitySolver(mesh_mean, bathy, xs, ef, flow_model, opts, 'NoExpand', V, flow_regs); 
[flow_solv.rotation] = deal(xs.angle);
flow = flow_solv.get_solution();

%% Empirical model: Taylor
flow_model = TaylorVelocityModel;
% flow_model.constituents = constituents;
flow_model.s_order = deal([1 1 1]);
flow_model.n_order = deal([1 1 1]);
flow_model.sigma_order = deal([1 1 1]);

%% Solver options and regularization

% %% Add regularization terms
flow_regs = regularization.Velocity.get_all_regs(mesh_mean, bathy, xs, flow_model, 'NoExpand', V);

flow_regs(1).weight = [10 1000];
flow_regs(2).weight = [10 1000];
flow_regs(3).weight = [1 100];
flow_regs(4).weight = [1 100];
flow_regs(5).weight = [10 1000];
%% Solve

flow_solv = LocationBasedVelocitySolver(mesh_mean, bathy, xs, ef, flow_model, opts, 'NoExpand', V, flow_regs); 
[flow_solv.rotation] = deal(xs.angle);
flow = flow_solv.get_solution();
%%
for cc = 1:numel(V)
    for reg = 1:2
     mesh_mean(cc).plot(flow(cc).pars(:,42, reg),...
          'AspectRatio',5, 'FixAspectRatio', true);
%     title(flow_model(cc).all_names(flow_model(cc).find_par(order = 0, component = 'u', variable = 'sig', constituent = 'M0')))
    colorbar;
    end
end
%%

for cc = 1:numel(V)
    for reg = 1:2
    figure;
     mesh_mean(cc).plot(flow(cc).pars(:,flow_model(cc).find_par(order = 0, component = 'u', variable = 'sig', constituent = 'M0'), reg),...
          'AspectRatio',5, 'FixAspectRatio', true);
    title(flow_model(cc).all_names(flow_model(cc).find_par(order = 0, component = 'u', variable = 'sig', constituent = 'M0')))
    colorbar;
%      nexttile;
     mesh_mean(cc).plot(flow(cc).pars(:,flow_model(cc).find_par(order = 0, component = 'u', variable = 'sig', constituent = 'M2', sincos = "cos"), reg),...
          'AspectRatio',5, 'FixAspectRatio', true);
    title(flow_model(cc).all_names(flow_model(cc).find_par(order = 0, component = 'u', variable = 'sig', constituent = 'M2', sincos = "cos")))

    colorbar;
%      nexttile;
     mesh_mean(cc).plot(flow(cc).pars(:,flow_model(cc).find_par(order = 0, component = 'u', variable = 'sig', constituent = 'M2', sincos = "sin"), reg),...
          'AspectRatio',5, 'FixAspectRatio', true);
         title(flow_model(cc).all_names(flow_model(cc).find_par(order = 0, component = 'u', variable = 'sig', constituent = 'M2', sincos = "sin")))

    colorbar;
    mesh_mean(cc).plot(flow(cc).pars(:,flow_model(cc).find_par(order = 1, component = 'u', variable = 'sig', constituent = 'M0'), 2),...
          'AspectRatio',5, 'FixAspectRatio', true);
        title(flow_model(cc).all_names(flow_model(cc).find_par(order = 1, component = 'u', variable = 'sig', constituent = 'M0')))

    colorbar;
    end
end

%% Backscatter
bran = 1; % I have 3 datasets, all contained in one VMADCP object.
bac_model(bran) = TidalScalarModel; % Type your own model here, Shiyu. In my case I use a tidalmodel, just to test.
[bac_model(bran).constituents] = deal(constituents);

bac_solv = BackscatterSolver(mesh_mean(bran), bathy(bran), xs(bran), ef(bran), bac_model(bran), opts, 'NoExpand', V);
bac = bac_solv.get_solution();
figure;
mesh_mean.plot(bac.pars(:,1)); colorbar;
% mesh_mean.plot(bac.pars(:,2)); colorbar;
% mesh_mean.plot(bac.pars(:,3)); colorbar;
% mesh_mean.plot(bac.pars(:,4)); colorbar;
% mesh_mean.plot(bac.pars(:,5)); colorbar;

%% Salinity

sal_model(bran) = TidalScalarModel;
[sal_model(bran).constituents] = deal(constituents);

% sal_regs = regularization.Scalar.get_all_regs(mesh_mean, bathy, xs, sal_model, 'NoExpand', V);

sal_solv = ExternalDataSolver(mesh = mesh_mean(bran), bathy = bathy(bran), xs=xs(bran), ef=ef(bran), data_model=sal_model(bran), opts=opts,...
    position = [ctd{bran}.X, ctd{bran}.Y, ctd{bran}.Z],...
    time = datetime(ctd{bran}.T2, 'ConvertFrom','datenum'),...
    data = ctd{bran}.S,...
    xform = 1,...
    water_level_object = water_level);

disp("made solver")

sal = sal_solv.get_solution();

%         position(:,3) double = [] % x,y,z position of the data
%         time(:,1) datetime = [] % time the data were measured
%         data(:,1) double = [] % data values
%         xform(:,:) double = [] % transformation matrix to transform the data
%         water_level_object(1,1) WaterLevel % water levels


% gof = MP.get_residuals();
% 
% MP.plot_residuals({0,1,5})
% 
% MP.plot_solution({data_model.names{1}{1:5}}, 1:3, 'v', 0, 'w', 0, 'ArrowScaling', [30, 20], 'ArrowTransform','symlog', 'ArrowParam', [.001, .1]);
% 
% MP.plot_solution({data_model.names{1}{6},...
%     data_model.names{2}{11},...
%     data_model.names{3}{16}}, 1:3, 'v', 0, 'w', 0, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);

%MP.plot_residuals();
% %% Visualize mesh
% m = makefigure(20,3);
% subplot(1,2,1)
% MP.plot_mesh('LineStyle', '-', 'FixAspectRatio', false, 'vertical', 'z');
% title('$z$ - coordinates', 'Interpreter','latex','FontSize', 12);
% xlabel([])
% 
% subplot(1,2,2)
% MP.plot_mesh('LineStyle', '-', 'FixAspectRatio', false, 'vertical', 'sig');
% xlabel([])
% %gca.TickLabelInterpreter = 'latex';
% 
% title('$\sigma$ - coordinates', 'Interpreter','latex', 'FontSize', 12);
% axs=axes(m,'visible','off');
% %m.Title.Visible='on';
% axs.XLabel.Visible='on';
% %m.YLabel.Visible='on';
% %ylabel(m,'yourYLabel');
% xlabel(axs,'y [m]', 'Interpreter','latex');
% %title(m,'yourTitle');
% %% Post Processing
% % inspect all
% 
% 
% 
% MP.plot_solution({data_model.names{1}{6:10}}, 1:3, 'v', 1, 'w', 1, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);
% 
% MP.plot_solution({data_model.names{1}{11:15}}, 1:3, 'v', 1, 'w', 1, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);
% 
% MP.plot_solution({data_model.names{1}{16:20}}, 1:3, 'v', 1, 'w', 1, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);
% 
% % Gradient terms
% 
% MP.plot_solution({data_model.names{1}{6},...
%     data_model.names{2}{11},...
%     data_model.names{3}{16}, ...
%     data_model.names{1}{7},...
%     data_model.names{2}{12},...
%     data_model.names{3}{17}}, 1:3, 'v', 0, 'w',0, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);
% 
% 
% MP.plot_solution({data_model.names{1}{6},...
%     data_model.names{2}{11},...
%     data_model.names{3}{16}}, 1:3, 'v', 0, 'w', 0, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);
% 
% 
% MP.plot_solution({data_model.names{1}{7},...
%     data_model.names{2}{12},...
%     data_model.names{3}{17}}, 1:3, 'v', 0, 'w', 0, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);
% 
% 
% 
% 
% %% Sensitivity
% %MP.cross_validation();
% 
% CV = MP.cross_validation_analysis();
% 
% MP.plot_cross_validation_analysis(CV, 1);
% 
% SA = MP.local_sensitivity_analysis();
