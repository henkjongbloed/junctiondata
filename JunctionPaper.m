%% ADCP data: PAPER 3
% Created by Henk Jongbloed

clearvars % Clear all variables from the workspace
%close all % Close all figures
figsave = 0;

RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'WP2\TwoJunctions\code'])

ds = 1;
tak = 1;

DS = {'OMHA14','NMOMNW15'}; DS = DS{ds};
BN = {{' Hartel Canal'; ' Old Meuse South'; ' Old Meuse North'},...
    {' New Meuse'; ' Old Meuse'; ' Rotterdam Waterway'}}; BN = BN{ds};


%Import preprocessed semi-raw data (ADCP / CTD / H structs)

[adcp, ctd, h] = import_data(RF, DS);

%% Make VMADCP

V = rdi.VMADCP(adcp{tak});
if strcmp(DS, 'NMOMNW15')
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


load(strcat("xs_", DS)); % define the cross-section
load(strcat("ef_", DS));
xs = xs(tak);
ef = ef(tak); % -> define further fitering such that nensembles equals number of ef ensembles

V.filters = Filter;
bathy = BathymetryScatteredPoints(ef, V);
bathy.interpolator.span = 0.01;

mesh_maker = SigmaZetaMeshFromVMADCP(ef, xs, bathy, 'NoExpand', V);

mesh_mean = mesh_maker.get_mesh(resn = 50, resz = 20); % get mesh at mean water level
% mesh_mean.plot()

%% SolverOptions
opts = SolverOptions(extrapolate_vert = 0);
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

flow_regs(1).weight = [1000];
flow_regs(2).weight = [1000];
flow_regs(3).weight = [100];
flow_regs(4).weight = [100];
flow_regs(5).weight = [1000];

%% Solve and plot

flow_solv = LocationBasedVelocitySolver(mesh_mean, bathy, xs, ef, flow_model, opts, 'NoExpand', V, flow_regs); 
flow_solv.rotation = xs.angle;
flow = flow_solv.get_solution();

flow.plot_solution({flow_model.all_names{[1]}}, representation = "Aphi", sol_idx = 1) % u,v,w
%%
nqt = 100; nqn = 50; nqs = 20;
t = linspace(0, 1, nqt); % in days
n = linspace(min(mesh_mean.n_left)+.5, max(mesh_mean.n_right)-.5, nqn);
sig = linspace(0, 1, nqs);

[T, N, Sig] = ndgrid(t,n,sig); Tq = T(:); Nq = N(:); Sigq = Sig(:);

U = flow.evaluate(1, T = Tq, N = Nq, Sig = Sigq, extrapolate=true);

u = reshape(U(:,1), [nqt, nqn, nqs]);
v = reshape(U(:,2), [nqt, nqn, nqs]);
w = reshape(U(:,3), [nqt, nqn, nqs]);

% Decompose

X.t = t; X.y = n; X.sig = sig;
X.T = T; X.Y = N; X.Sig = Sig;


[Hf, Wlf, Zbf] = get_H(X, mesh_mean, bathy, water_level, 1); % Own function to derive H(t,y) from X and adcptools variables.

H = repmat(Hf, [1,1,numel(X.sig)]);
Zb = repmat(Zbf, [1,1,numel(X.sig)]); 
X.Z = Zbf + X.Sig.*H;

flow.decomp = Decomposition(X=X, H=H);

[DF, AF, DFf, AFf] = flow.decomp.decompose_function(u);

flow.decomp.plot_components(AF)

Dsum = sum(cat(4, DFf{:}), 4); %time x lateral x vertical
%% plot
%Dsum = H;
% fi = figure;
% filename = 'testAnimated.gif';
% ncolor = 100;
% amax = max(abs([min(Dsum, [], "all"), max(Dsum, [], "all")]));
% levels = linspace(-amax, amax, ncolor);
% [~,ha]=contourf(squeeze(X.Y(1,:,:))', squeeze(X.Z(1,:,:))', squeeze(Dsum(1,:,:))' , levels, "LineColor",'none');
% colorbar;
% title("Hartelkanaal: Flow")
% colormap(gca, flipud(brewermap(ncolor, 'RdBu')));
% clim([-amax, amax])
% ylim([min(X.Z, [], 'all'), max(X.Z, [], 'all')])
% set(gca, 'XDir','reverse') % Very important
% for tim = 1:1:nqt
%     frame = getframe(fi);
%     im = frame2im(frame);
%     ha.YData = squeeze(X.Z(tim,:,:))';
%     ha.ZData = squeeze(Dsum(tim,:,:))';
%     if tim == 1
%         [imind,cm] = rgb2ind(im,256);
%         imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
%     else
%         imind = rgb2ind(im, cm);
%         imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
%     end
%     pause(.2)
% end
% 

[SF,C] = flow.decomp.decompose_product(u,v);
figure
bar(SF)
% figure
flow.decomp.plotC(C)
%bar(SF)



%% Salinity

if strcmp(DS,'NMOMNW15')    % Change ADCP coordinates rather than CTD coordinates!!
      [ctd{tak}.pos(:,1), ctd{tak}.pos(:,2)]      = wgs2utm(ctd{tak}.Y,ctd{tak}.X, 31, 'N'); %Convert LatLon to UTM %Correct
      time = rescale(ctd{tak}.T, datenum(V.time(1)), datenum(V.time(end)));
      ctd{tak}.t = datetime(time, 'ConvertFrom', 'datenum'); % Rescale time such that measurements fall within experiment time
      ctd{tak}.pos(:,3) = water_level.get_water_level(ctd{tak}.t) - ctd{tak}.Z;
%     [ctdH, S.i] = getHour(ctd.T2, hour/2);% + datenum('9-Aug-2014 00:00:00'); % Correct
%     S.t = rescale(ctdH, min(U.adcpH), max(U.adcpH)); % Everything works perfectly
%     % Space
      xnew = center_cloud([ctd{tak}.pos(:,1), ctd{tak}.pos(:,2)], xs.origin'); % Because of offset in CTD coordinates relative to ADCP - quick fix
      ctd{tak}.pos(:,1) = xnew(:,1); ctd{tak}.pos(:,2) = xnew(:,2); % Only needed for 2015 data
      
elseif strcmp(DS,'OMHA14')
    [ctd{tak}.lon, ctd{tak}.lat, ~]      = rd2wgs(ctd{tak}.X, ctd{tak}.Y); %Convert RD to LatLon %Correct
    [ctd{tak}.pos(:,1), ctd{tak}.pos(:,2)]      = wgs2utm(ctd{tak}.lat, ctd{tak}.lon, 31, 'N'); %Convert LatLon to UTM %Correct
    time = rescale(ctd{tak}.T, datenum(V.time(1)), datenum(V.time(end)));
    ctd{tak}.t = datetime(time, 'ConvertFrom', 'datenum'); % Rescale time such that measurements fall within experiment time
    ctd{tak}.pos(:,3) = -ctd{tak}.Z;% Upward facing z axis - max(-ctd{tak}.Z - water_level.get_water_level(ctd{tak}.t));
    [ctd{tak}.s, ctd{tak}.n] = xs(tak).xy2sn(ctd{tak}.pos(:,1), ctd{tak}.pos(:,2));
    %     ctddt = 6000;
%     [ctdH, S.i] = getHour(ctd.T, ctddt); % dt = 6000 -> 1 sample per 6s?
     %if strcmp(U.BN, ' Old Meuse North')
%         [S.t, S.i] = manualOMN(ctd.T, U.adcpH);
%     end
end

Xs.t = time;
Xs.x = ctd{tak}.pos(:,1);
Xs.y = ctd{tak}.pos(:,2);
[Hctd, Wlctd, Zbctd]  = get_H(Xs, mesh_mean, bathy, water_level, 0);
hctd = diag(Hctd);
zbctd = diag(Zbctd);
%zpos = ctd{tak}.pos(:,3) + min(ctd{tak}.pos(:,3));
sigctd = (ctd{tak}.pos(:,3)-zbctd)./hctd; % Rescale and scale back to z
sig_corr = rescale(sigctd, 0, 1);
plotctd(ctd{tak}, mesh_mean, water_level)

ctd{tak}.pos(:,3) = sigctd.*hctd + zbctd;

figure;
mesh_mean.plot3()
hold on
scatter3(ctd{tak}.pos(:,1), ctd{tak}.pos(:,2), ctd{tak}.pos(:,3))


opts = SolverOptions(extrapolate_vert = 1);
sal_model = TidalScalarModel;
sal_model.constituents = constituents;

sal_regs = regularization.Scalar.get_all_regs(mesh_mean, bathy, xs, sal_model, 'NoExpand', V);
sal_regs(1).weight = 1;
sal_regs(2).weight = 1;
sal_solv = ExternalDataSolver(mesh=mesh_mean, bathy=bathy, xs=xs, model=sal_model, opts=opts, regularization=sal_regs, ...
    position = ctd{tak}.pos,...
    time = ctd{tak}.t,...
    data = ctd{tak}.S,...
    xform = ones(size(ctd{tak}.S)),...
    water_level_object = water_level);

% disp("made solver")

sal = sal_solv.get_solution();
figure;
mesh_mean.plot(sal.pars(:,1)); colorbar;
mesh_mean.plot(sal.pars(:,2)); colorbar;
mesh_mean.plot(sal.pars(:,3)); colorbar;
mesh_mean.plot(sal.pars(:,4)); colorbar;
mesh_mean.plot(sal.pars(:,5)); colorbar;

S = sal.evaluate(1, T = Tq, N = Nq, Sig = Sigq, extrapolate=true);
S(S<0)=0;
s = reshape(S(:,1), [nqt, nqn, nqs]);
sal.decomp = flow.decomp;
[DFs, AFs, DFfs, AFfs] = sal.decomp.decompose_function(s);
sal.decomp.plot_components(AFs)
[SF,C] = sal.decomp.decompose_product(u,s);

figure
bar(SF)
% figure
sal.decomp.plotC(C)
bar(SF)


Dsums = sum(cat(4, DFfs{:}), 4); %time x lateral x vertical
%% plot
% Salinity: positive

%figure;
fi = figure;
filename = 'testAnimatedSal.gif';
ncolor = 100;
%Dsums = H;
amax = max(abs([min(Dsums, [], "all"), max(Dsums, [], "all")]));
levels = linspace(0, amax, ncolor);
[~,ha]=contourf(squeeze(X.Y(1,:,:))', squeeze(X.Z(1,:,:))', squeeze(Dsums(1,:,:))' , levels, "LineColor",'none');
colorbar;
title("Hartelkanaal: Salinity")
cm = flipud(brewermap(ncolor, 'RdBu'));
colormap(gca, cm(size(cm,1)/2:end,:));
clim([0, amax])
ylim([min(X.Z, [], 'all'), max(X.Z, [], 'all')])
set(gca, 'XDir','reverse') % Very important
for tim = 1:1:nqt
    ha.YData = squeeze(X.Z(tim,:,:))';
    ha.ZData = squeeze(Dsums(tim,:,:))';
    frame = getframe(fi);
    im = frame2im(frame);
    if tim == 1
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
    else
        imind = rgb2ind(im, cm);
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
    end
    pause(.1)
end

%% Plot
% flow.plot_solution(flow_model.all_names(1:5))

%% Plot

% for cc = 1:numel(V)
%     for reg = 1:2
%      mesh_mean(cc).plot(flow(cc).pars(:,1, reg),...
%           'AspectRatio',5, 'FixAspectRatio', true);
% %     title(flow_model(cc).all_names(flow_model(cc).find_par(order = 0, component = 'u', variable = 'sig', constituent = 'M0')))
%     colorbar;
%     end
% end




%% Backscatter
% tak = 1; % I have 3 datasets, all contained in one VMADCP object.
% bac_model(tak) = TidalScalarModel;
% [bac_model(tak).constituents] = deal(constituents);
% bac_regs = regularization.Scalar.get_all_regs(mesh_mean, bathy, xs, bac_model, 'NoExpand', V);
% 
% bac_regs(1).weight = 1;
% bac_regs(2).weight = 1;
% 
% bac_solv = BackscatterSolver(mesh_mean(tak), bathy(tak), xs(tak), ef(tak), bac_model(tak), opts, 'NoExpand', V, bac_regs);
% bac = bac_solv.get_solution();
% 
% figure;
% mesh_mean.plot(bac.pars(:,1)); colorbar;
% mesh_mean.plot(bac.pars(:,2)); colorbar;
% mesh_mean.plot(bac.pars(:,3)); colorbar;
% mesh_mean.plot(bac.pars(:,4)); colorbar;
% mesh_mean.plot(bac.pars(:,5)); colorbar;

%         position(:,3) double = [] % x,y,z position of the data
%         time(:,1) datetime = [] % time the data were measured
%         data(:,1) double = [] % data values
%         xform(:,:) double = [] % transformation matrix to transform the data
%         water_level_object(1,1) WaterLevel % water levels


