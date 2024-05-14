%% ADCP data: Method paper
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

BNs = {{'HK2'; 'OMS2'; 'OMN2'},...
    {'NM2'; 'OM2'; 'RW2'}}; BNs = BNs{ds};

%Import preprocessed semi-raw data (ADCP / CTD / H structs)

[adcp, ctd, h] = import_data(RF, DS);

%% Make VMADCP
for t = 1:numel(tak)
     V(t) = rdi.VMADCP(adcp{tak(t)}); 
     if strcmp(DS, 'NMOMNW15')
         V(t).horizontal_position_provider(1) = [];
     end
end

%% Initial inspection
constituents = {'M2', 'M4'};
% 
water_level = VaryingWaterLevel(datetime(datevec(h.t)), h.wl); % set water level to varying
water_level.model = TidalScalarModel(constituents = constituents);
water_level.model.scalar_name = 'eta'; % Scalar
water_level.get_parameters();
[V.water_level_object] = deal(water_level);



%% create bathymetry
%bathy = BathymetryScatteredPoints(V); % create a bathymetry (includes all data)

dat_selection = "saved";
switch dat_selection
    case "manual"
        [ef, xs] = cross_section_selector(V);
    case "default"
        xs = XSection(V); % define the cross-section
        ef = EnsembleFilter(V);
    case "saved" % Save ensemblefilter and cross-section once
        load(strcat("xs_", DS)); % define the cross-section
        load(strcat("ef_", DS));
        xs = xs(tak);
        ef = ef(tak); % -> define further fitering such that nensembles equals number of ef ensembles
    otherwise 
        xs = XSection(V); % define the cross-section
        ef = EnsembleFilter(V);
end

[V.filters] = deal(Filter);
bathy = BathymetryScatteredPoints(ef, V);
interpolators = [bathy.interpolator];
[interpolators.span] = deal(0.01);

mesh_maker = SigmaZetaMeshFromVMADCP(ef, xs, bathy, 'NoExpand', V);

mesh_mean = mesh_maker.get_mesh(resn = 50, resz = 20); % get mesh at mean water level
% mesh_mean.plot()


% %% Mesh plot
% % mesh_mean.plot('FixAspectRatio', false)
% 
% %% Visualize mesh
% % m = makefigure(10,3);
% % mesh_mean(1).plot('sigma', false, 'FixAspectRatio', false);
% % title('$z$ - coordinates', 'Interpreter','latex','FontSize', 12);
% % axis tight
% % ax = gca;
% % set(ax, 'XDir','reverse')
% % ax.TickLabelInterpreter = 'latex';
% % axs=axes(m,'visible','off');
% % %m.Title.Visible='on';
% % axs.XLabel.Visible='on';
% % axs.YLabel.Visible='on';
% % %ylabel(m,'yourYLabel');
% % xlabel(axs,'y [m]', 'Interpreter','latex');
% % ylabel(axs, '$z [m]$', 'Interpreter', 'latex')
% % 
% % 
% % m = makefigure(10,3);
% % mesh_mean(1).plot('sigma', true, 'FixAspectRatio', true);
% % title('$\sigma$ - coordinates', 'Interpreter','latex', 'FontSize', 12);
% % axis tight
% % ax = gca;
% % set(ax, 'XDir','reverse')
% % ax.TickLabelInterpreter = 'latex';
% % axs=axes(m,'visible','off');
% % %m.Title.Visible='on';
% % axs.XLabel.Visible='on';
% % axs.YLabel.Visible='on';
% % %ylabel(m,'yourYLabel');
% % xlabel(axs,'y [m]', 'Interpreter','latex');
% % ylabel(axs,'$\sigma$', 'Interpreter', 'latex')
% 
%% SolverOptions
opts = SolverOptions(extrapolate_vert = 0);
opts.force_zero = [1 1 1 1 1];


%% Empirical model: TaylorTidal
flow_model = TaylorTidalVelocityModel;
flow_model.constituents = constituents;
flow_model.s_order = deal([1 1 1]);
flow_model.n_order = deal([1 1 1]);
flow_model.sigma_order = deal([1 1 1]);

%% Solver options and regularization

% %% Add regularization terms
flow_regs = regularization.Velocity.get_all_regs(mesh_mean, bathy, xs, flow_model, opts, 'NoExpand', V);

flow_regs(1).weight = [1000];
flow_regs(2).weight = [1000];
flow_regs(3).weight = [100];
flow_regs(4).weight = [100];
flow_regs(5).weight = [1000];

% flow_regs(1).weight = [1000*3^2];
% flow_regs(2).weight = [1000*3^2];
% flow_regs(3).weight = [100*861^2];
% flow_regs(4).weight = [100*3^2];
% flow_regs(5).weight = [1000*.3^2];
%% Solve and plot

flow_solv = LocationBasedVelocitySolver(mesh_mean, bathy, xs, ef, flow_model, opts, 'NoExpand', V, flow_regs); 
[flow_solv.rotation] = deal(xs.angle);
flow = flow_solv.get_solution();

%%
figure;
mesh_mean.plot(flow.pars(:,1)); colorbar;
mesh_mean.plot(flow.pars(:,2)); colorbar;
mesh_mean.plot(flow.pars(:,3)); colorbar;
mesh_mean.plot(flow.pars(:,4)); colorbar;
mesh_mean.plot(flow.pars(:,5)); colorbar;

figure;
mesh_mean.plot(flow.pars(:,21)); colorbar;
mesh_mean.plot(flow.pars(:,22)); colorbar;
mesh_mean.plot(flow.pars(:,23)); colorbar;
mesh_mean.plot(flow.pars(:,24)); colorbar;
mesh_mean.plot(flow.pars(:,25)); colorbar;

figure;
mesh_mean.plot(flow.pars(:,41)); colorbar;
mesh_mean.plot(flow.pars(:,42)); colorbar;
mesh_mean.plot(flow.pars(:,43)); colorbar;
mesh_mean.plot(flow.pars(:,44)); colorbar;
mesh_mean.plot(flow.pars(:,45)); colorbar;


%%
nqt = 100; nqn = 50; nqs = 20;
t = linspace(0, 1, nqt); % in days
n = linspace(min(mesh_mean.n_left)+.5, max(mesh_mean.n_right)-.5, nqn);%*ones(size(t));
sig = linspace(0, 1, nqs);%*ones(size(t)); %TODO let sig go fron 0 to 1

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

%H = 6*ones(size(X.T));

flow.decomp = Decomposition(X=X, H=H);

[DF, AF, DFf, AFf] = flow.decomp.decompose_function(u);

flow.decomp.plot_components(AF)

%D.animate()

Dsum = sum(cat(4, DFf{:}), 4); %time x lateral x vertical
%% plot
%Dsum = H;
fi = figure;
filename = 'testAnimated.gif';
ncolor = 100;
amax = max(abs([min(Dsum, [], "all"), max(Dsum, [], "all")]));
levels = linspace(-amax, amax, ncolor);
[~,ha]=contourf(squeeze(X.Y(1,:,:))', squeeze(X.Z(1,:,:))', squeeze(Dsum(1,:,:))' , levels, "LineColor",'none');
colorbar;
title("Hartelkanaal: Flow")
colormap(gca, flipud(brewermap(ncolor, 'RdBu')));
clim([-amax, amax])
ylim([min(X.Z, [], 'all'), max(X.Z, [], 'all')])
set(gca, 'XDir','reverse') % Very important
for tim = 1:1:nqt
    frame = getframe(fi);
    im = frame2im(frame);
    ha.YData = squeeze(X.Z(tim,:,:))';
    ha.ZData = squeeze(Dsum(tim,:,:))';
    if tim == 1
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
    else
        imind = rgb2ind(im, cm);
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
    end
    pause(.2)
end


% [SF,C] = D.decompose_product(u,v);
% figure
% bar(SF)
% % figure
% D.plotC(C)
% %bar(SF)



%% Salinity
t=1;
if strcmp(DS,'NMOMNW15')    % Change ADCP coordinates rather than CTD coordinates!!
      [ctd{tak(t)}.pos(:,1), ctd{tak(t)}.pos(:,2)]      = wgs2utm(ctd{tak(t)}.Y,ctd{tak(t)}.X, 31, 'N'); %Convert LatLon to UTM %Correct
      time = rescale(ctd{tak(t)}.T, datenum(V.time(1)), datenum(V.time(end)));
      ctd{tak(t)}.t = datetime(time, 'ConvertFrom', 'datenum'); % Rescale time such that measurements fall within experiment time
      ctd{tak(t)}.pos(:,3) = water_level.get_water_level(ctd{tak(t)}.t) - ctd{tak(t)}.Z;
%     [ctdH, S.i] = getHour(ctd.T2, hour/2);% + datenum('9-Aug-2014 00:00:00'); % Correct
%     S.t = rescale(ctdH, min(U.adcpH), max(U.adcpH)); % Everything works perfectly
%     % Space
      xnew = center_cloud([ctd{tak(t)}.pos(:,1), ctd{tak(t)}.pos(:,2)], xs.origin'); % Because of offset in CTD coordinates relative to ADCP - quick fix
      ctd{tak(t)}.pos(:,1) = xnew(:,1); ctd{tak(t)}.pos(:,2) = xnew(:,2); % Only needed for 2015 data
      
elseif strcmp(DS,'OMHA14')
    [ctd{tak(t)}.lon, ctd{tak(t)}.lat, ~]      = rd2wgs(ctd{tak(t)}.X, ctd{tak(t)}.Y); %Convert RD to LatLon %Correct
    [ctd{tak(t)}.pos(:,1), ctd{tak(t)}.pos(:,2)]      = wgs2utm(ctd{tak(t)}.lat, ctd{tak(t)}.lon, 31, 'N'); %Convert LatLon to UTM %Correct
    time = rescale(ctd{tak(t)}.T, datenum(V.time(1)), datenum(V.time(end)));
    ctd{tak(t)}.t = datetime(time, 'ConvertFrom', 'datenum'); % Rescale time such that measurements fall within experiment time
    ctd{tak(t)}.pos(:,3) = -ctd{tak(t)}.Z;% Upward facing z axis - max(-ctd{tak(t)}.Z - water_level.get_water_level(ctd{tak(t)}.t));
    [ctd{tak(t)}.s, ctd{tak(t)}.n] = xs(tak(t)).xy2sn(ctd{tak(t)}.pos(:,1), ctd{tak(t)}.pos(:,2));
    %     ctddt = 6000;
%     [ctdH, S.i] = getHour(ctd.T, ctddt); % dt = 6000 -> 1 sample per 6s?
     %if strcmp(U.BN, ' Old Meuse North')
%         [S.t, S.i] = manualOMN(ctd.T, U.adcpH);
%     end
end

%figure;




Xs.t = time;
Xs.x = ctd{tak(t)}.pos(:,1);
Xs.y = ctd{tak(t)}.pos(:,2);
[Hctd, Wlctd, Zbctd]  = get_H(Xs, mesh_mean, bathy, water_level, 0);
hctd = diag(Hctd);
zbctd = diag(Zbctd);
%zpos = ctd{tak(t)}.pos(:,3) + min(ctd{tak(t)}.pos(:,3));
sigctd = (ctd{tak(t)}.pos(:,3)-zbctd)./hctd; % Rescale and scale back to z
sig_corr = rescale(sigctd, 0, 1);
plotctd(ctd{tak(t)}, mesh_mean, water_level)

ctd{tak(t)}.pos(:,3) = sigctd.*hctd + zbctd;

figure;
mesh_mean.plot3()
hold on
scatter3(ctd{tak(t)}.pos(:,1), ctd{tak(t)}.pos(:,2), ctd{tak(t)}.pos(:,3))
% S.eta = U.eta(S.t); % Waterlevel per CTD sample

% S.Z          =  - ctd.Z + S.eta;

% [S.s, S.n] = U.xs.xy2sn(S.X, S.Y);

% S.zb = n2zb(S.n, U.mesh_mean);
% S.sig = z2sig(S.Z, S.zb, S.eta); % IS THIS ALLOWED IN SIGMA COORDINATES?? - TAKE CARE
% S.Z = sig2z(S.sig, S.zb, S.eta);

% [S.Xp, S.Yp] = U.xs.sn2xy(0*S.s, S.n); % Project on cross section by setting s = 0. [Xp, Yp] = projected X,Y

%      max(adcpT) - min(adcpT)
%      max(ctd.T) - min(ctd.T)
%     % Continue here: Make S.t corresponding to 
%     S.t = ctd.T;% + datenum('12-Aug-2014 00:00:00');
    %U.mesh_mean.plot()


opts = SolverOptions(extrapolate_vert = 0);
sal_model = TidalScalarModel;
[sal_model.constituents] = deal(constituents);

sal_regs = regularization.Scalar.get_all_regs(mesh_mean, bathy, xs, sal_model, 'NoExpand', V);
sal_regs(1).weight = 1;
sal_regs(2).weight = 1;
sal_solv = ExternalDataSolver('NoExpand', regularization = sal_regs, mesh=mesh_mean, bathy=bathy, xs=xs, ef=ef, data_model=sal_model, opts=opts,...
    position = ctd{tak(t)}.pos,...
    time = ctd{tak(t)}.t,...
    data = ctd{tak(t)}.S,...
    xform = ones(size(ctd{tak(t)}.S)),...
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
%bar(SF)


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

% figure
% plot(u(:,1, 1));
% figure
% contourf(squeeze(u(:,:,5))');
% xlabel('t'); ylabel('y')
% colorbar


% figure
% contourf(squeeze(u(:,10,:))');
% xlabel('t'); ylabel('sig')
% colorbar
%%
% figure
% for tim = 1:nqt
%     subplot(3,1,1)
%     contourf(squeeze(u(tim,:,:))');
%     colorbar
%     clim([-max(abs(u), [],  "all"), max(abs(u),[], "all")])
% 
%     subplot(3,1,2)
%     contourf(squeeze(v(tim,:,:))');
%     colorbar
%     clim([-max(abs(v),[],  "all"), max(abs(v),[], "all")])
% 
%     subplot(3,1,3)
%     contourf(squeeze(w(tim,:,:))');
%     colorbar
%     clim([-max(abs(w), [],  "all"), max(abs(w),[], "all")])
%     pause(.2)
%     
% end
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
% tak(t) = 1; % I have 3 datasets, all contained in one VMADCP object.
% bac_model(tak(t)) = TidalScalarModel;
% [bac_model(tak(t)).constituents] = deal(constituents);
% bac_regs = regularization.Scalar.get_all_regs(mesh_mean, bathy, xs, bac_model, 'NoExpand', V);
% 
% bac_regs(1).weight = 1;
% bac_regs(2).weight = 1;
% 
% bac_solv = BackscatterSolver(mesh_mean(tak(t)), bathy(tak(t)), xs(tak(t)), ef(tak(t)), bac_model(tak(t)), opts, 'NoExpand', V, bac_regs);
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
% %% Visualize mesh
% m = makefigure(20,3);
% subplot(1,2,1)
% MP.plot_mesh('LineStyle', '-', 'FixAspectRatio', false, 'vertical', 'z');
% title('$z$ - coordinates', 'Interpreter','latex','FontSize', 12);
% xlabel([])
% 
% subplot(1,2,2)
% MP.plot_mesh('LineStyle', '-', 'FixAspectRatio', false, 'vertical', 'sig');
% xlabel([])
% %gca.TickLabelInterpreter = 'latex';
% 
% title('$\sigma$ - coordinates', 'Interpreter','latex', 'FontSize', 12);
% axs=axes(m,'visible','off');
% %m.Title.Visible='on';
% axs.XLabel.Visible='on';
% %m.YLabel.Visible='on';
% %ylabel(m,'yourYLabel');
% xlabel(axs,'y [m]', 'Interpreter','latex');
% %title(m,'yourTitle');
% %% Post Processing
% % inspect all
% 
% 
% 
% MP.plot_solution({data_model.names{1}{6:10}}, 1:3, 'v', 1, 'w', 1, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);
% 
% MP.plot_solution({data_model.names{1}{11:15}}, 1:3, 'v', 1, 'w', 1, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);
% 
% MP.plot_solution({data_model.names{1}{16:20}}, 1:3, 'v', 1, 'w', 1, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);
% 
% % Gradient terms
% 
% MP.plot_solution({data_model.names{1}{6},...
%     data_model.names{2}{11},...
%     data_model.names{3}{16}, ...
%     data_model.names{1}{7},...
%     data_model.names{2}{12},...
%     data_model.names{3}{17}}, 1:3, 'v', 0, 'w',0, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);
% 
% 
% MP.plot_solution({data_model.names{1}{6},...
%     data_model.names{2}{11},...
%     data_model.names{3}{16}}, 1:3, 'v', 0, 'w', 0, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);
% 
% 
% MP.plot_solution({data_model.names{1}{7},...
%     data_model.names{2}{12},...
%     data_model.names{3}{17}}, 1:3, 'v', 0, 'w', 0, 'ArrowScaling', [30, 2], 'ArrowTransform','symlog', 'ArrowParam', [.001, .001]);
% 
% 
% 
% 
% %% Sensitivity
% %MP.cross_validation();
% 
% CV = MP.cross_validation_analysis();
% 
% MP.plot_cross_validation_analysis(CV, 1);
% 
% SA = MP.local_sensitivity_analysis();


% for cc = 1:numel(V)
%     for reg = 1:2
%     figure;
%      mesh_mean(cc).plot(flow(cc).pars(:,flow_model(cc).find_par(order = 0, component = 'u', variable = 'sig', constituent = 'M0'), reg),...
%           'AspectRatio',5, 'FixAspectRatio', true);
%     title(flow_model(cc).all_names(flow_model(cc).find_par(order = 0, component = 'u', variable = 'sig', constituent = 'M0')))
%     colorbar;
% %      nexttile;
%      mesh_mean(cc).plot(flow(cc).pars(:,flow_model(cc).find_par(order = 0, component = 'u', variable = 'sig', constituent = 'M2', sincos = "cos"), reg),...
%           'AspectRatio',5, 'FixAspectRatio', true);
%     title(flow_model(cc).all_names(flow_model(cc).find_par(order = 0, component = 'u', variable = 'sig', constituent = 'M2', sincos = "cos")))
% 
%     colorbar;
% %      nexttile;
%      mesh_mean(cc).plot(flow(cc).pars(:,flow_model(cc).find_par(order = 0, component = 'u', variable = 'sig', constituent = 'M2', sincos = "sin"), reg),...
%           'AspectRatio',5, 'FixAspectRatio', true);
%          title(flow_model(cc).all_names(flow_model(cc).find_par(order = 0, component = 'u', variable = 'sig', constituent = 'M2', sincos = "sin")))
% 
%     colorbar;
%     mesh_mean(cc).plot(flow(cc).pars(:,flow_model(cc).find_par(order = 1, component = 'u', variable = 'sig', constituent = 'M0'), 2),...
%           'AspectRatio',5, 'FixAspectRatio', true);
%         title(flow_model(cc).all_names(flow_model(cc).find_par(order = 1, component = 'u', variable = 'sig', constituent = 'M0')))
% 
%     colorbar;
%     end
% end
