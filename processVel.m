function [V, U] = processVel(ADCP, h, BN)

V = VMADCP(ADCP); % create VMADCP object
V.horizontal_position_provider(1) = [];
V.water_level_object = VaryingWaterLevel(datetime(datevec(h.t)), h.wl); % set water level to varying
% plot(h.t, h.wl)
%disp(V.water_level)
%% create bathymetry
B = BathymetryScatteredPoints(V); % create a bathymetry (includes all data)
%B.plot
%B.plot_residuals % have a good look at the residuals of the bathymetry.

xs = XSection(V); % define the cross-section
if any(strcmp({' New Meuse'; ' Rotterdam Waterway';' Hartel Canal'}, BN))
    xs.revert();
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
mesh_maker.deltaz = (Hw+1)/7;

mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level
% mesh_mean.plot()
% hold on
% plot(mesh_mean.n_middle(mesh_mean.col_to_cell), mesh_mean.z_center, '*')
% hold off
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
T = LocationBasedVelocitySolver(V, xs, ef, mesh_mean, B);


T.velocity_model = TidalVelocityModel();
TM2 = 12.42; %hours

T.velocity_model.constituentsU = [TM2; TM2/2];
T.velocity_model.constituentsV = [TM2; TM2/2];
T.velocity_model.constituentsW = [TM2; TM2/2];
[U.pars, U.cov_pars, U.n_bvels] = T.get_parameters;

[U.pars{1,1}, U.cov_pars{1,1}] = T.rotate_to_xs_pars(U.pars{1,1}, U.cov_pars{1,1});

% surf(abs(U.pars{1,1} - U.pars_xy{1,1}))

%% Here: Rotate to direction of u0, v0 such that uo // xs_orth. Make new cross section.

% n = length(T.velocity_model.constituentsU); %Dirty trick
% idx = [1 2:2:2*n];
% % U.pars{1,1}(:,idx) = - U.pars{1,1}(:,idx);
% U.pars{1,1}(:,1) = - U.pars{1,1}(:,1);

[U.tid_pars, ~] = T.velocity_model.get_tidal_pars(U.pars{1,1}, U.cov_pars{1,1});

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
end


U.adcpT = datenum(V.raw.timeV ) ; % in days
hour = datenum('12-Aug-2014 01:00:00') - datenum('12-Aug-2014 00:00:00');

[U.adcpH, U.adcpI] = getHour(U.adcpT, hour/2); % ADCP hours

U.salmap = brewermap(15, 'YlOrBr'); 
U.velmap = brewermap(20, 'RdBu');
U.phimap = [U.velmap ;flipud(U.velmap)];
U.Bw = Bw;

end


