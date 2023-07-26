%% ADCP data: Method paper
% Created by Henk Jongbloed

clearvars % Clear all variables from the workspace
close all % Close all figures
figsave = 0;

RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'WP2\TwoJunctions\code'])

ds = 1;
tak = 2;

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
mesh_maker.deltan = (Bw + 1)/25;%(Bw + 1)/35; -> 35 ill-conditioned -> reg trivially works better.
mesh_maker.deltaz = (Hw+1)/13;%(Hw+1)/13; -> 13 ill-conditioned -> reg trivially works better.
mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level


%% Mesh plot
mesh_mean.plot('FixAspectRatio', false)

%% Empirical model

data_model = TaylorTidalModel(constituents = constituents, s_order = [1,1,1], n_order = [1,1,1], sigma_order = [1,1,1]);

opts = SolverOptions(reg_iter = [3,3,3,3,3], max_reg_pars = [1000, 1000, 400, 400, 1000]);
opts.algorithm = "pcg";
opts.pcg_iter = 4000;
opts.cv_iter = 1;
opts.pcg_tol = 1e-8;
opts.reg_pars = {0*[1,1,1,1,1], 3*[10,10,1,1,10], 25*[10,10,10,10,10]};
%% Regularization
R = Regularization(bathy = bathy, xs = xs, mesh = mesh_mean, model = data_model);
R.assemble_matrices()

data_model.rotation = xs.angle();
S = LocationBasedVelocitySolver(V, mesh_mean, bathy, xs, ef, data_model, R, opts);
MP = S.get_parameters();
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
%% Visualize mesh
m = makefigure(20,3);
subplot(1,2,1)
MP.plot_mesh('LineStyle', '-', 'FixAspectRatio', false, 'vertical', 'z');
title('$z$ - coordinates', 'Interpreter','latex','FontSize', 12);
xlabel([])

subplot(1,2,2)
MP.plot_mesh('LineStyle', '-', 'FixAspectRatio', false, 'vertical', 'sig');
xlabel([])
%gca.TickLabelInterpreter = 'latex';

title('$\sigma$ - coordinates', 'Interpreter','latex', 'FontSize', 12);
axs=axes(m,'visible','off');
%m.Title.Visible='on';
axs.XLabel.Visible='on';
%m.YLabel.Visible='on';
%ylabel(m,'yourYLabel');
xlabel(axs,'y [m]', 'Interpreter','latex');
%title(m,'yourTitle');
%% Post Processing
% inspect all



MP.plot_solution({data_model.names{1}{6:10}}, 1:3, 'v', 1, 'w', 1, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);

MP.plot_solution({data_model.names{1}{11:15}}, 1:3, 'v', 1, 'w', 1, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);

MP.plot_solution({data_model.names{1}{16:20}}, 1:3, 'v', 1, 'w', 1, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);

% Gradient terms

MP.plot_solution({data_model.names{1}{6},...
    data_model.names{2}{11},...
    data_model.names{3}{16}, ...
    data_model.names{1}{7},...
    data_model.names{2}{12},...
    data_model.names{3}{17}}, 1:3, 'v', 0, 'w',0, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);


MP.plot_solution({data_model.names{1}{6},...
    data_model.names{2}{11},...
    data_model.names{3}{16}}, 1:3, 'v', 0, 'w', 0, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);


MP.plot_solution({data_model.names{1}{7},...
    data_model.names{2}{12},...
    data_model.names{3}{17}}, 1:3, 'v', 0, 'w', 0, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);




%% Sensitivity
%MP.cross_validation();

CV = MP.cross_validation_analysis();

MP.plot_cross_validation_analysis(CV, 1);

SA = MP.local_sensitivity_analysis();
