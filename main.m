%% ADCP data: Method paper
% Created by Henk Jongbloed

clearvars % Clear all variables from the workspace
close all % Close all figures
figsave = 0;

RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'WP2\TwoJunctions\code'])

ds = 1;
tak = 1:3

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
% 
% V.plot_orientations;
% figure
% V.plot_track
% figure
% V.plot_bed_position
% figure
% V.plot_track_velocity
% figure
% V.plot_velocity
% figure
% V.plot_backscatter
% create VMADCP object
% for cc = 1:numel(V)
%     V(cc).horizontal_position_provider(1) = [];
%     mtc(cc) = MagneticDeviationTwoCycle;
%     mtc(cc).plot_tracks(V(cc))
%     mtc(cc).estimate_deviation(V(cc))
%     mtc(cc).estimate_deviation(V(cc),true)
% end
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
bathy = BathymetryScatteredPoints(ef, 'NoExpand', V);
interpolators = [bathy.interpolator];
[interpolators.span] = deal(0.01);

mesh_maker = SigmaZetaMeshFromVMADCP(ef, xs, bathy,'NoExpand', V);
[mesh_maker.deltan] = deal(8);
[mesh_maker.deltaz] = deal(.8);
mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level


% Bw = mesh_mean.nw(2) - mesh_mean.nw(1); % Little fishy but it works
% Hw = abs(min(mesh_mean.zb_all));
% mesh_maker.deltan = (Bw + 1)/25;%(Bw + 1)/35; -> 35 ill-conditioned -> reg trivially works better.
% mesh_maker.deltaz = (Hw+1)/13;%(Hw+1)/13; -> 13 ill-conditioned -> reg trivially works better.
% mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level
% 

%% Mesh plot
% mesh_mean.plot('FixAspectRatio', false)

%% Empirical model
%Trying different models

flow_model = TaylorTidalVelocityModel;


flow_model.constituents = deal(constituents);
flow_model.s_order = deal([1 1 1]);
flow_model.n_order = deal([1 1 1]);
flow_model.sigma_order = deal([1 1 1]);

opts = SolverOptions(reg_iter = [3,3,3,3,3], max_reg_pars = [1000, 1000, 400, 400, 1000]);
opts.algorithm = "pcg";
opts.pcg_iter = 4000;
opts.cv_iter = 1;
opts.pcg_tol = 1e-8;

opts = SolverOptions;
flow_regs = regularization.Velocity.get_all_regs(mesh_mean, bathy, xs, flow_model, 'NoExpand', V);
% regs = regularization.Regularization(mesh_mean, bathy, xs, data_model, 'NoExpand', V);


for cc = 1:numel(V) % discuss with bart
    disp([flow_regs{cc}.weight])
    [flow_regs{cc}.weight] = deal(0);
    disp([flow_regs{cc}.weight])
   % disp(flow_regs{cc}(1, 3).weight)
end

flow_solv = LocationBasedVelocitySolver(mesh_mean, bathy, xs, ef, flow_model, opts, flow_regs, 'NoExpand', V); 
[flow_solv.rotation] = deal(xs.angle);
flow = flow_solv.get_solution();

figure
for cc = 1:numel(V)
        mesh_mean(cc).plot(flow(cc).pars(:,1),'AspectRatio',5, 'FixAspectRatio', true);
%     mesh_mean(cc).plot(mp(cc).pars(:,data_model(cc).find_par(1, {'u','v','w'}, {'sig'}, {'M0'})),...
%         'AspectRatio',5, 'FixAspectRatio', true);
    colorbar;
end

%% Salinity
bran = 1;
sal_model(bran) = TidalScalarModel;
[sal_model(bran).constituents] = deal(constituents);

sal_regs = regularization.Scalar.get_all_regs(mesh_mean, bathy, xs, sal_model, 'NoExpand', V);

sal_solv = ExternalDataSolver(mesh_mean(bran), bathy(bran), xs(bran), ef(bran), sal_model(bran), opts,...
    position = [ctd{bran}.X, ctd{bran}.Y, ctd{bran}.Z],...
    time = datetime(ctd{bran}.T, 'ConvertFrom','datenum'),...
    data = ctd{bran}.S,...
    xform = 1,...
    water_level_object = water_level);

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
