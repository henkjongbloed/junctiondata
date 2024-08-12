function [flow, salt] = solve_flow_salinity(dname, bi, adcp, ctd, h, lf, ls, ures, sres)
%% Make VMADCP

V = rdi.VMADCP(adcp{bi});
if strcmp(dname, 'NMOMNW15')
    V.horizontal_position_provider(1) = [];
end


%% Initial inspection
constituents = {'M2', 'M4'};
%

tlim = [V.time(1); V.time(end)];
dt = 1/24; % two hours
[tc, idxtc] = clip_time(datetime(datevec(h.t)), tlim, dt);


water_level = VaryingWaterLevel(tc, h.wl(idxtc)); % set water level to varying
water_level.model = TidalScalarModel(constituents = constituents);
water_level.model.scalar_name = 'eta'; % Scalar
water_level.get_parameters();
V.water_level_object = water_level;

%% create bathymetry
%bathy = BathymetryScatteredPoints(V); % create a bathymetry (includes all data)


xs = load(strcat("xs_", dname)).xs(bi); % define the cross-section
ef = load(strcat("ef_", dname)).ef(bi);

xs.revert; % Revert all XSections

V.filters = Filter;
bathy = BathymetryScatteredPoints(ef, V);
bathy.interpolator.span = 0.01;

mesh_maker = SigmaZetaMeshFromVMADCP(ef, xs, bathy, 'NoExpand', V);

mesh_mean = mesh_maker.get_mesh(resn = ures(1), resz = ures(2)); % get mesh at mean water level
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


mesh_mean_s = mesh_maker.get_mesh(resn = sres(1), resz = sres(2)); % get mesh at mean water level


ctd = prep_ctd(dname, ctd, bi, V, water_level, xs, mesh_mean_s, bathy);

sopts = SolverOptions(extrapolate_vert = 1, lat_weight_factor = 10); % ramp up requirement of lateral homogeneity
sal_model = TidalScalarModel;
sal_model.constituents = constituents;
sal_model.scalar_name = "salt";
sal_regs = regularization.Scalar.get_all_regs(mesh_mean_s, bathy, xs, sal_model, 'NoExpand', V);
sal_regs(1).weight = ls(:,1)';
sal_regs(2).weight = ls(:,2)';
sal_solv = ExternalDataSolver(mesh=mesh_mean_s, bathy=bathy, xs=xs, model=sal_model, opts=sopts, regularization=sal_regs, ...
    position = ctd{bi}.pos,...
    time = ctd{bi}.t,...
    data = ctd{bi}.S,...
    xform = ones(size(ctd{bi}.S)),...
    water_level_object = water_level);

% disp("made solver")

salt = sal_solv.get_solution();


end