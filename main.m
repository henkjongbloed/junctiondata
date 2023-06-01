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

V = rdi.VMADCP(adcp{tak}); % create VMADCP object
% V.horizontal_position_provider(1) = [];
mtc= MagneticDeviationTwoCycle;
%mtc.plot_tracks(V)
%mtc.estimate_deviation(V)
%mtc.estimate_deviation(V,true)

constituents = {'M2', 'M4'};
% 
water_level = VaryingWaterLevel(datetime(datevec(h.t)), h.wl); % set water level to varying
water_level.model = TidalModel(constituents = constituents);
water_level.model.components = {'eta'}; % Scalar
water_level.get_parameters();

V.water_level_object = water_level;
%% create bathymetry
bathy = BathymetryScatteredPoints(V); % create a bathymetry (includes all data)

dat_selection = "saved";
switch dat_selection
    case "manual"
        [ef, xs] = cross_section_selector(V);
    case "default"
        xs = XSection(V); % define the cross-section
        ef = EnsembleFilter(V);
    case "saved" % Save ensemblefilter and cross-section once
        load("xs_HA14"); % define the cross-section
        load("ef_HA14");
    otherwise 
        xs = XSection(V); % define the cross-section
        ef = EnsembleFilter(V);
end

% if any(strcmp({' New Meuse'; ' Rotterdam Waterway';' Hartel Canal'}, BN))
%     xs.revert();
% end

% U.BN = BN;
% 
% if any(strcmp({' New Meuse'; ' Old Meuse'; ' Rotterdam Waterway'}, BN))
%     U.DS = 2015;
% else
%     U.DS = 2014;
% end


%% create mesh

mesh_maker = SigmaZetaMeshFromVMADCP(V, bathy, xs); % mesh generator
mesh_maker.deltan = 100;
mesh_maker.deltaz = 2;
mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level
Bw = mesh_mean.nw(2) - mesh_mean.nw(1); % Little fishy but it works
Hw = abs(min(mesh_mean.zb_all));
mesh_maker.deltan = (Bw + 1)/25;%(Bw + 1)/35;
mesh_maker.deltaz = (Hw+1)/11;%(Hw+1)/13;
mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level

data_model = TaylorTidalModel(constituents = constituents, s_order = [1,1,1], n_order = [1,1,1], sigma_order = [1,1,1]);

opts = SolverOptions();
opts.algorithm = "pcg";
opts.pcg_iter = 4000;
opts.cv_iter = 1;
opts.pcg_tol = 1e-8;
opts.reg_pars = {0*[1,1,1,1,1], 3*[10,10,1,1,10], 100*[10,10,1,1,10]};

R = Regularization(bathy = bathy, xs = xs, mesh = mesh_mean, model = data_model);
R.assemble_matrices()

data_model.rotation = xs.angle();
S = LocationBasedVelocitySolver(V, mesh_mean, bathy, xs, ef, data_model, R, opts);
MP = S.get_parameters();

% MP.cross_validate();

figure;
MP.plot_solution({data_model.names{1}{1:5}}, 1:3, 'v', 0, 'w', 0, 'ArrowScaling', [10, 10], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);

MP.cross_validation();

MP.cross_validation_analysis()
