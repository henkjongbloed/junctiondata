restoredefaultpath
clearvars
% close all

addpath ~/src/adcptools/ %https://github.com/bartverm/adcptools   -> point this to the right path

%% load data
dat_file = 'raw_data.mat';
if exist(dat_file,"file")
    load(dat_file)
else
    raw_data=rdi.readDeployment('trans','raw_data');
    save(dat_file,'raw_data')
end

vmadcp=rdi.VMADCP(raw_data);

def_file = "transect_defs.mat";
if exist(def_file,"file")
    load(def_file)
else
    [ef, xs]=cross_section_selector(vmadcp);
    save(def_file,"xs","ef")
end
%  % select cross-sections. now its loaded in next line

% xs.revert;

%% Correct heading
Mtc = MagneticDeviationTwoCycle;
Mtc.plot_tracks(vmadcp)
Mtc.estimate_deviation(vmadcp,true);
vmadcp.heading_provider(2).magnetic_deviation_model = Mtc;
Mtc.plot_tracks(vmadcp)


%% process repeat transects
vmadcp.filters = Filter;
bathy = BathymetryScatteredPoints(vmadcp,ef);
interpolators = [bathy.interpolator];
[interpolators.span] = deal(0.01);

mesh_maker=SigmaZetaMeshFromVMADCP(vmadcp,ef,xs,bathy);
mesh = mesh_maker.get_mesh();

data_model = TaylorVelocityModel( ...
    s_order = [1 0 1],...
    n_order = [0 1 1],...
    sigma_order = [1 1 1]);

regs = regularization.Velocity.get_all_regs(mesh, bathy, xs, data_model, vmadcp);
opts = SolverOptions;
solv = LocationBasedVelocitySolver(mesh, bathy, xs, ef, data_model, regs, opts, vmadcp); 
% [solv(1).regularization.weight] = deal(100, 100, 5, 5, 100);
[solv.rotation] = deal(xs.angle);

mp = solv.get_parameters();

%% plot
figure
for cc = 1:7
    mesh(cc).plot(mp(cc).pars(:,data_model.find_par(0,{'u','v','w'})),...
        'AspectRatio',5, 'FixAspectRatio', true)
    colorbar
end

