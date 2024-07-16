function [flow, sal] = solve_flow_salinity(adcp, ctd, tak_idx, dataset_idx, h, lf, ls, res)
%% Make VMADCP

V = rdi.VMADCP(adcp{tak_idx});
if strcmp(dataset_idx, 'NMOMNW15')
    V.horizontal_position_provider(1) = [];
end


%% Initial inspection
constituents = {'M2', 'M4'};
% 
water_level = VaryingWaterLevel(datetime(datevec(h.t)), h.wl); % set water level to varying
water_level.model = TidalScalarModel(constituents = constituents);
water_level.model.scalar_name = 'eta'; % Scalar
water_level.get_parameters();
V.water_level_object = water_level;

%% create bathymetry
%bathy = BathymetryScatteredPoints(V); % create a bathymetry (includes all data)


xs = load(strcat("xs_", dataset_idx)).xs(tak_idx); % define the cross-section
ef = load(strcat("ef_", dataset_idx)).ef(tak_idx);

V.filters = Filter;
bathy = BathymetryScatteredPoints(ef, V);
bathy.interpolator.span = 0.01;

mesh_maker = SigmaZetaMeshFromVMADCP(ef, xs, bathy, 'NoExpand', V);

mesh_mean = mesh_maker.get_mesh(resn = res(1), resz = res(2)); % get mesh at mean water level
% mesh_mean.plot()

%% SolverOptions
opts = SolverOptions(extrapolate_vert = 0, lat_weight_factor = 10);
opts.force_zero = [1 1 1 1 1];


%% Empirical model: TaylorTidal
flow_model = TaylorTidalVelocityModel;
flow_model.constituents = constituents;
flow_model.s_order = [1 1 1];
flow_model.n_order = [1 1 1];
flow_model.sigma_order = [1 1 1];

%% Solver options and regularization

% %% Add regularization terms
flow_regs = regularization.Velocity.get_all_regs(mesh_mean, bathy, xs, flow_model, opts, 'NoExpand', V);

flow_regs(1).weight = lf(:,1)';
flow_regs(2).weight = lf(:,2)';
flow_regs(3).weight = lf(:,3)';
flow_regs(4).weight = lf(:,4)';
flow_regs(5).weight = lf(:,5)';

%% Solve for the flow

flow_solv = LocationBasedVelocitySolver(mesh_mean, bathy, xs, ef, flow_model, opts, 'NoExpand', V, flow_regs); 
flow_solv.rotation = xs.angle;
flow = flow_solv.get_solution();

%% Solve for the salinity
ctd = prep_ctd(dataset_idx, ctd, tak_idx, V, water_level, xs, mesh_mean, bathy);

sopts = SolverOptions(extrapolate_vert = 1, lat_weight_factor = 10);
sal_model = TidalScalarModel;
sal_model.constituents = constituents;
sal_model.scalar_name = "salt";
sal_regs = regularization.Scalar.get_all_regs(mesh_mean, bathy, xs, sal_model, 'NoExpand', V);
sal_regs(1).weight = ls(:,1)';
sal_regs(2).weight = ls(:,2)';
sal_solv = ExternalDataSolver(mesh=mesh_mean, bathy=bathy, xs=xs, model=sal_model, opts=sopts, regularization=sal_regs, ...
    position = ctd{tak_idx}.pos,...
    time = ctd{tak_idx}.t,...
    data = ctd{tak_idx}.S,...
    xform = ones(size(ctd{tak_idx}.S)),...
    water_level_object = water_level);

% disp("made solver")

sal = sal_solv.get_solution();


end