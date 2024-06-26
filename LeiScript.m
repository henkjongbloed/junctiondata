%% Script Lei
load('C:\Users\jongb013\Documents\PHD\2-Programming\WP2\TwoJunctions\data\LeiData\raw_CYP3.mat')
ADCP = raw_dat;
load('C:\Users\jongb013\Documents\PHD\2-Programming\WP2\TwoJunctions\data\LeiData\B-CYP3.mat')

load('C:\Users\jongb013\Documents\PHD\2-Programming\WP2\TwoJunctions\data\LeiData\M.mat')
load('C:\Users\jongb013\Documents\PHD\2-Programming\WP2\TwoJunctions\data\LeiData\MESH-GEN.mat')
%load('C:\Users\jongb013\Documents\PHD\2-Programming\WP2\TwoJunctions\data\LeiData\VMADCP.mat')

V = VMADCP(ADCP); % create VMADCP object
V.horizontal_position_provider(1) = [];
V.water_level_object = ConstantWaterLevel;%  VaryingWaterLevel(datetime(datevec(h.t)), h.wl); % set water level to varying
%V.water_level_object.constituents = {};
%V.water_level_object.get_parameters();

% V.water_level_object.get_omega();


% plot(h.t, h.wl)
%disp(V.water_level)
%% create bathymetry
%B = BathymetryScatteredPoints(V); % create a bathymetry (includes all data)
% B.plot

%B.plot_residuals % have a good look at the residuals of the bathymetry.
xs = XSection(V); % define the cross-section

mesh_mean = M;

%% create mesh

mesh_mean.plot()
hold on
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
if strcmp(model, 'tidal')
    T.velocity_model = TidalVelocityModel();

    T.velocity_model.constituentsU = {'M2', 'M4'};
    T.velocity_model.constituentsV = {'M2', 'M4'};
    T.velocity_model.constituentsW = {'M2', 'M4'};
    [U.pars, U.cov_pars, U.n_bvels] = T.get_parameters;

    % [U.pars{1,1}, U.cov_pars{1,1}] = T.rotate_to_xs_pars(U.pars{1,1}, U.cov_pars{1,1});

    % surf(abs(U.pars{1,1} - U.pars_xy{1,1}))

    %% Here: Rotate to direction of u0, v0 such that uo // xs_orth. Make new cross section.

    % n = length(T.velocity_model.constituentsU); %Dirty trick
    % idx = [1 2:2:2*n];
    % % U.pars{1,1}(:,idx) = - U.pars{1,1}(:,idx);
    % U.pars{1,1}(:,1) = - U.pars{1,1}(:,1);

    [U.tid_pars, U.cov_tid_pars] = T.velocity_model.get_tidal_pars(U.pars{1,1}, U.cov_pars{1,1});
elseif strcmp(model, 'steadyTaylor')
    T.velocity_model = TaylorVelocityModel();
    T.velocity_model.s_order = [1,1,0];
    T.velocity_model.n_order = [1,1,0];
    T.velocity_model.sigma_order = [1,1,1];
    [U.pars, U.cov_pars, U.n_bvels] = T.get_parameters;

elseif strcmp(model, 'taylorTidal')
    T.velocity_model = TaylorTidalVelocityModel();

    % Initiate model
    T.velocity_model.s_order = [1,1,1];
    T.velocity_model.n_order = [1,1,1];
    T.velocity_model.sigma_order = [1,1,1];
    T.velocity_model.z_order = [0,0,0];
    T.velocity_model.constituentsU = {};
    T.velocity_model.constituentsV = {};
    T.velocity_model.constituentsW = {};
    %T.velocity_model.constituentsU = {'M2', 'M4'};
    %T.velocity_model.constituentsV = {'M2', 'M4'};
    %T.velocity_model.constituentsW = {'M2', 'M4'};
    % Reg_Pars
    opts.reg_vary = 'coupled'; % coupled, full
    opts.reg_iter = [12,nan,12,nan, nan]; %if reg_vary neq none
    opts.res_near_zero = .05; % Smaller: more resulution near zero.
    opts.max_reg_pars = 1e3; %max reg pars 
    ymax = symlog(opts.max_reg_pars, opts.res_near_zero);

    opts.reg_pars = {symexp(linspace(0,ymax,opts.reg_iter(1)), opts.res_near_zero)',...
        symexp(linspace(0,ymax,opts.reg_iter(2)), opts.res_near_zero)',...
        symexp(linspace(0,ymax,opts.reg_iter(3)), opts.res_near_zero)',...
        symexp(linspace(0,ymax,opts.reg_iter(4)), opts.res_near_zero)',...
        symexp(linspace(0,ymax,opts.reg_iter(5)), opts.res_near_zero)'};

    % In the above, keep 0 as fixed value!
    opts.reg_pars0 = {100, 100, 5, 5, 100}; % Default reg parameters. Can also be multiple for purposes of comparing.
    opts.reg_pars1 = {[0, 100, 1000], [0, 100, 1000], [0, 20, 50], [0, 20, 50], [0, 100, 1000]}; % Default reg parameters. Can also be multiple for purposes of comparing.


    opts.set_zero = [0,0,0,0,0]; % Set to one if you want to set lambda_i to zero.

    % CV_Pars
    opts.cv_mode = 'none'; %none, random, omit_cells, omit_time
    opts.omit_cells = 42:59; % in case of omit_cells
    opts.training_perc = .90; % in case of random
    opts.cv_iter = 5; % Cross-validation iterations (using different set partitions) in case of random
    opts.use_p0 = 1;
    % Solve
    [tempp, tempcp, tempnbv] = T.get_parameters_reg(opts);
    [tempp0, tempcp0, tempnbv0] = T.get_parameters;
    [U.pars, U.cov_pars, U.n_bvels] = deal(tempp{1,1}, tempcp{1,1}, tempnbv{1,1});
    [U.pars0, U.cov_pars0, U.n_bvels0] = deal(tempp0{1,1}, tempcp0{1,1}, tempnbv0{1,1});
end

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

