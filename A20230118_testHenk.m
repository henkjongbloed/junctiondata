% Test script to read and process data with method Henk

fclose all; clearvars; close all; clc
restoredefaultpath

%% Add paths
RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\';
RF2 = 'WP2\TwoJunctions\';
DATADIR = [RF, RF2, 'data\20221121_iris'];
c_path = [DATADIR, filesep, '\CROSS']; l_path = [DATADIR, filesep, '\LONG'];
addpath(c_path);
addpath(l_path);

addpath(strcat(RF,'Tools\adcptoolsGit')); %path to ADCPTools
addpath(strcat(RF,'Tools\loess-master')); %path to loess-master
addpath(strcat(RF,'Tools\BrewerMap-master')); %path to BrewerMap (see github)

%% Settings

date = '20210308';
% freq = 1200;
freq = 600;
readnew = 1;

%% Read data

% Read water levels
%wl = readwldata(['Data\waterlevels_', date, '.csv'], {'Maassluis'}); wl = wl{1,1};
load('C:\Users\jongb013\Documents\PHD\2-Programming\WP2\TwoJunctions\data\20221121_iris/wl_20210308_voorHenk.mat')
% Determine location of adcp data
DATLOC = l_path;

% Start reading ADCP data
if readnew      % read adcp data
    allcrs = dir(fullfile([DATLOC, '\A*r.000']));
    for nt = 1:length(allcrs)
        PATH = fullfile(DATLOC);
        adcp{1,nt} = readDeployment(allcrs(nt).name, PATH);
    end
    nloc = size(adcp,1);

    % Add spatial information
    for nt = 1:length(allcrs)
        x{1,nt}=[];y{1,nt}=[];
        if ~isempty(adcp{1,nt})
            for n=1:size(adcp{1,nt}.timeV,1)
                adcp{1,nt}.time1(n) = datenum(adcp{1,nt}.timeV(n,:));
            end
            [tx, ty] = utmADCP(adcp{1,nt});
            x{1,nt} = tx;
            y{1,nt} = ty;
            clear tx ty
        end
    end

    % adcp_full = readDeployment('NW', [cd,DATLOC]);
    adcp_full = readDeployment('A', DATLOC );

    %     % Save
    %     save([date, '_adcp', num2str(freq), '_ALL.mat'], 'adcp', 'ldb*', 'freq', 'date', 'nloc', 'wl', 'x', 'y');
else
    %     load([date, '_adcp', num2str(freq), '_CROSS.mat'])
    %load([date, '_adcp', num2str(freq), '_LONG.mat'])
    load('C:\Users\jongb013\Documents\PHD\2-Programming\WP2\TwoJunctions\data\20221121_iris\adcp600_FULL_20210308.mat')
end

%cmap = polarmap(64);    % Define colormap for certain figures
f=25;                   % Scale factor for quivers

%% Convert to type ADCP
rtime = datetime(datestr(wl.Time));

%for n = 1:length(adcp_full)
V = rdi.VMADCP(adcp_full);
%     VM{n}.type = ADCP_Type.RioGrande;

wl_int = interp1(rtime, wl.Val, V.time);
water_level = VaryingWaterLevel(V.time, wl_int);
V.water_level_object = water_level;
clear water_level wl_int
%end

% A_full = ADCP(adcp_full);
% V = VMADCP(adcp_full);

%%
V.water_level_object.constituents = {'M2', 'M4'};
V.water_level_object.get_parameters();
V.water_level_object.get_omega();


% plot(h.t, h.wl)
%disp(V.water_level)
%% create bathymetry
B = BathymetryScatteredPoints(V); % create a bathymetry (includes all data)
% B.plot

%B.plot_residuals % have a good look at the residuals of the bathymetry.
xs = XSection(V); % define the cross-section
% if any(strcmp({' New Meuse'; ' Rotterdam Waterway';' Hartel Canal'}, BN))
xs.revert();
% end
hold on
xs.plot()
hold off

% U.BN = BN;
%
% if any(strcmp({' New Meuse'; ' Old Meuse'; ' Rotterdam Waterway'}, BN))
%     U.DS = 2015;
% else
%     U.DS = 2014;
% end


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

mesh_maker.deltan = (Bw + 1)/11;
mesh_maker.deltaz = (Hw+1)/5;

mesh_mean = mesh_maker.get_mesh(); % get mesh at mean water level
% mesh_mean.get_neighbors(cell_idx)
% figure;
% plotMap(mesh_mean, U.DS)
figure;
mesh_mean.plot()
%hold on
% plot(mesh_mean.n_middle(mesh_mean.col_to_cell), mesh_mean.z_center, '*')

for c = 1:mesh_mean.ncells
    [nb, dom] = mesh_mean.get_neighbors(c);

    txt = {[c, dom]};
    text(mesh_mean.n_middle(mesh_mean.col_to_cell(c)), mesh_mean.z_center(c), txt)
end
hold off
%mesh_mean.plot()


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
    T.velocity_model.constituentsU = {'M2', 'M4'};
    T.velocity_model.constituentsV = {'M2', 'M4'};
    T.velocity_model.constituentsW = {'M2', 'M4'};

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
%wl = @(t) interp1(h.t, h.wl, t, 'linear');

U.B = B; U.xs = xs; U.mesh_mean = mesh_mean;
% U.mesh_mean = mesh_mean;
U.T = T;
% U.vel = vel;
% U.cov_vel = cov_vel;
U.eta = wl;

% U.BN = BN;
%
% if any(strcmp({' New Meuse'; ' Old Meuse'; ' Rotterdam Waterway'}, BN))
%     U.DS = 2015;
% else
%     U.DS = 2014;
% end


U.adcpT = datenum(V.raw.timeV ) ; % in days
hour = datenum('12-Aug-2014 01:00:00') - datenum('12-Aug-2014 00:00:00');



[U.adcpH, U.adcpI] = getHour(U.adcpT, hour/2); % ADCP hours

U.salmap = brewermap(15, 'YlOrBr');
U.velmap = brewermap(20, 'RdBu');
U.phimap = [U.velmap ;flipud(U.velmap)];
U.Bw = Bw;

names = U.T.velocity_model.names;
%             p0 = pab2pAphi(dat.p0, U{1}.T.velocity_model.names);
%         makefigure(50, 22)
%          plotp(dat.p1, names, names,  10, U{1})
[p1(:,1), tid_names] = pab2pAphi(dat.p1(:,1), U.T.velocity_model.names);
p1(:,2) = pab2pAphi(dat.p1(:,2), U.T.velocity_model.names);
p1(:,3) = pab2pAphi(dat.p1(:,3), U.T.velocity_model.names);

figure;
plotp_compare({p1(:,1), p1(:,2), p1(:,3)}, tid_names([21,22,23,24,25,36,6,31,56]), tid_names,  U)

figure;
plotp_compare({dat.p1(:,1), dat.p1(:,2), dat.p1(:,3)},names([1,2,3,4,5,16,6,31,56]), names,  U)
figure;
plotp(p1(:,2), tid_names, tid_names,  6, U)
figure;
plotp(dat.p1(:,2), names, names,  6, U)


