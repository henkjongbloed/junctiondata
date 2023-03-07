function [V, U] = processVel(ADCP, h, BN)

V = VMADCP(ADCP); % create VMADCP object
V.horizontal_position_provider(1) = [];
V.water_level_object = VaryingWaterLevel(datetime(datevec(h.t)), h.wl); % set water level to varying
V.water_level_object.constituents = {'M2', 'M4'};
V.water_level_object.get_parameters();

V.water_level_object.get_omega();


% plot(h.t, h.wl)
%disp(V.water_level)
%% create bathymetry
B = BathymetryScatteredPoints(V); % create a bathymetry (includes all data)
%load('C:\Users\jongb013\Documents\PHD\2-Programming\WP2\TwoJunctions\data\LeiData\B-CYP3.mat')
% B.plot

%B.plot_residuals % have a good look at the residuals of the bathymetry.
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
mesh_maker = SigmaZetaMeshFromVMADCP(V, B, xs); % mesh generator

mesh_maker.deltan = 100;
mesh_maker.deltaz = 2;

% mesh_maker.deltan = 12;
% mesh_maker.deltaz = 1.5;


mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level
% mesh_mean.plot()

% if strcmp(BN, ' OMN')
%     mesh_maker.deltan = 3;
% end
%
% if strcmp(BN, ' OM')
%     mesh_maker.deltan = 5;
% end

Bw = mesh_mean.nw(2) - mesh_mean.nw(1); % Little fishy but it works
Hw = abs(min(mesh_mean.zb_all));

%mesh_maker.deltan = (Bw + 1)/12;
%mesh_maker.deltaz = 2*Hw/8;

mesh_maker.deltan = (Bw + 1)/25;
mesh_maker.deltaz = (Hw+1)/13;

mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level
% mesh_mean.get_neighbors(cell_idx)
% figure;
% plotMap(mesh_mean, U.DS)
% figure;
% mesh_mean.plot()
% hold on
% plot(mesh_mean.n_middle(mesh_mean.col_to_cell), mesh_mean.z_center, '*')

for c = 1:mesh_mean.ncells
    [nb, dom] = mesh_mean.get_neighbors(c);

    txt = {[c, dom]};
    text(mesh_mean.n_middle(mesh_mean.col_to_cell(c)), mesh_mean.z_center(c), txt)
end
hold off
% mesh_mean.plot()


%% split repeat_transects
ef = EnsembleFilter(V);
% % now expand to water level of other repeat transects
% wl = repmat(V.water_level, numel(ef),1); % get water levels and repeat for each transect
% wl(vertcat(ef.bad_ensembles)) = nan; % set water_levels not in transect to nan
% wl = mean(wl,2,"omitnan"); % compute average water level for each transect
% mesh = mesh_mean.mesh_at_water_level(wl); % compute expanded meshes

%%
% T = LocationBasedVelocitySolver(V, xs, ef, mesh_mean, B);
TM2 = 12.42; %hours

T = LocationBasedVelocitySolver(V, xs, ef, mesh_mean);
model = 'taylorTidal';
T.velocity_model = TaylorTidalVelocityModel();

% Initiate model
T.velocity_model.s_order = [1,1,1];
T.velocity_model.n_order = [1,1,1];
T.velocity_model.sigma_order = [1,1,1];
T.velocity_model.z_order = [0,0,0];
%T.velocity_model.constituentsU = {};
%T.velocity_model.constituentsV = {};
%T.velocity_model.constituentsW = {};
T.velocity_model.constituentsU = {'M2', 'M4'};
T.velocity_model.constituentsV = {'M2', 'M4'};
T.velocity_model.constituentsW = {'M2', 'M4'};
% Reg_Pars

% Solve
[tempp, tempcp, tempnbv] = T.get_parameters_reg(opts);
[tempp0, tempcp0, tempnbv0] = T.get_parameters;
[U.pars, U.cov_pars, U.n_bvels] = deal(tempp{1,1}, tempcp{1,1}, tempnbv{1,1});
[U.pars0, U.cov_pars0, U.n_bvels0] = deal(tempp0{1,1}, tempcp0{1,1}, tempnbv0{1,1});


% U.pars = U.tid_pars;
% U.cov_pars = U.cov_tid_pars;
wl = @(t) interp1(h.t, h.wl, t, 'linear');

U.B = B; U.xs = xs; U.mesh_mean = mesh_mean;
% U.mesh_mean = mesh_mean;
U.T = T;
% U.vel = vel;
% U.cov_vel = cov_vel;
U.eta = wl;

U.BN = BN;

if any(strcmp({' New Meuse'; ' Old Meuse'; ' Rotterdam Waterway'}, BN))
    U.DS = 2015;
else
    U.DS = 2014;
end


U.adcpT = datenum(V.raw.timeV ) ; % in days
hour = datenum('12-Aug-2014 01:00:00') - datenum('12-Aug-2014 00:00:00');

[U.adcpH, U.adcpI] = getHour(U.adcpT, hour/2); % ADCP hours

U.salmap = brewermap(15, 'YlOrBr');
U.velmap = brewermap(20, 'RdBu');
U.phimap = [U.velmap ;flipud(U.velmap)];
U.Bw = Bw;

end


