function [V, U] = processVel(ADCP, h, BN)

V = rdi.VMADCP(ADCP); % create VMADCP object
V.horizontal_position_provider(1) = [];

water_level = VaryingWaterLevel(datetime(datevec(h.t)), h.wl); % set water level to varying
water_level.model = TidalModel(constituents = {'M2', 'M4'});
water_level.model.components = {'eta'}; % Scalar
water_level.get_parameters();

V.water_level_object = water_level;
%% create bathymetry
bathy = BathymetryScatteredPoints(V); % create a bathymetry (includes all data)

xs = XSection(V); % define the cross-section
if any(strcmp({' New Meuse'; ' Rotterdam Waterway';' Hartel Canal'}, BN))
    xs.revert();
end
hold on
xs.plot()
hold off

U.BN = BN;

if any(strcmp({' New Meuse'; ' Old Meuse'; ' Rotterdam Waterway'}, BN))
    U.DS = 2015;
else
    U.DS = 2014;
end


%% create mesh
mesh_maker = SigmaZetaMeshFromVMADCP(V, bathy, xs); % mesh generator
mesh_maker.deltan = 100;
mesh_maker.deltaz = 2;
mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level
Bw = mesh_mean.nw(2) - mesh_mean.nw(1); % Little fishy but it works
Hw = abs(min(mesh_mean.zb_all));
mesh_maker.deltan = (Bw + 1)/25;
mesh_maker.deltaz = (Hw+1)/13;
mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level
%% split repeat_transects
% ef = EnsembleFilter(V);

%%

%tid = TidalModel(constituents = {'M2', 'M4'});
%tay = TaylorModel(s_order = [1,1,1], n_order = [1,1,1], sigma_order = [1,1,1]);

%data_model = TaylorTidalModel(taylor = tay, tidal = tid, rotation = pi/2);

data_model = TaylorTidalModel(constituents = {'M2', 'M4'}, s_order = [1,1,1], n_order = [1,1,1], sigma_order = [1,1,1]);



%R = Regularization(bathy, xs, mesh_mean, data_model);

R = [ContinuityRegularization(bathy, xs, mesh_mean, data_model, water_level);...
    CoherenceRegularization(bathy, xs, mesh_mean, data_model, water_level);...
    KinematicRegularization(bathy, xs, mesh_mean, data_model, water_level)];

S = Solver(mesh_mean, bathy, xs, R, V.water_level_object);




data_model.get_model()
data_model.rotate_matrix()
data_model.rotate_model()
S.get_parameters();



end