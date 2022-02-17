function [V, U] = processVel(ADCP, h, BN)

V = VMADCP(ADCP); % create VMADCP object
V.horizontal_position_provider(1) = [];
V.water_level_object = VaryingWaterLevel(datetime(datevec(h.time)),h.level); % set water level to varying
%disp(V.water_level)
%% create bathymetry
B = BathymetryScatteredPoints(V); % create a bathymetry (includes all data)
%B.plot
%B.plot_residuals % have a good look at the residuals of the bathymetry.

xs = XSection(V); % define the cross-section

%% create mesh
mesh_maker = SigmaZetaMeshFromVMADCP(V, B, xs); % mesh generator
mesh_maker.deltan = 100;
mesh_maker.deltaz = 2;
mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level

Bw = mesh_mean.nw(2) - mesh_mean.nw(1); % Little fishy but it works
Hw = abs(min(mesh_mean.zb_all));

%mesh_maker.deltan = (Bw + 1)/12;
%mesh_maker.deltaz = 2*Hw/8;

mesh_maker.deltan = (Bw + 1)/12;
mesh_maker.deltaz = 2*Hw;

mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level

%% split repeat_transects
% ef = manual_repeat_transects(V); %made function myself %Change this!!
ef = EnsembleFilter(V);
% % now expand to water level of other repeat transects
% wl = repmat(V.water_level, numel(ef),1); % get water levels and repeat for each transect
% wl(vertcat(ef.bad_ensembles)) = nan; % set water_levels not in transect to nan
% wl = mean(wl,2,"omitnan"); % compute average water level for each transect
% mesh = mesh_mean.mesh_at_water_level(wl); % compute expanded meshes

%%
% T = LocationBasedVelocitySolver(V, xs, ef, mesh_mean, B);
T = TimeBasedVelocitySolver(V, xs, ef, mesh_mean, B);

%L = LocationBasedVelocitySolver(V, xs, ef, mesh, B);

T.velocity_model = TidalVelocityModel();
TM2 = 12.42; %hours

T.velocity_model.constituentsU = [TM2];
T.velocity_model.constituentsV = [];%[TM2, TM2/2];
T.velocity_model.constituentsW = [];%[TM2, TM2/2];
[U.pars, U.cov_pars, U.n_bvels] = T.get_parameters;
[U.tid_pars, ~] = T.velocity_model.get_tidal_pars(U.pars{1,1}, U.cov_pars{1,1});

[vel, cov_vel] = T.velocity_model.get_velocity(U.pars{1,1}, U.cov_pars{1,1});

% nOW - GET_TIDAL

% [vel, cov_vel] = T.get_velocity; % 

%vel2 = L.get_velocity;




U.B = B; U.xs = xs; U.mesh_mean = mesh_mean; 
U.mesh = mesh_mean; 
U.T = T; 
U.vel = vel;
U.cov_vel = cov_vel;
%U.eta = wl;

U.BN = BN; 


end


