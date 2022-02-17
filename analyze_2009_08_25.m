restoredefaultpath
clearvars
addPaths
%addpath ~/src/adcptools-proctrans/
%addpath ~/src/loess/
% raw_data = readDeployment('trans','C:\Users\jongb013\OneDrive - WageningenUR\PhDHenkJongbloed\2. Programming\DataTools\Tools\example_adcp_processing\muaramuntai25Aug');
% vmadcp = VMADCP(raw_data);
% [ef, xs]=cross_section_selector(vmadcp);
%%
load data_2009_08_25.mat

%%
cs = 2;
for cs=1:numel(ef)
    bathy(cs) = BathymetryScatteredPoints(vmadcp,ef(cs));
    bathy(cs).interpolator.span = 0.01;
    mesh_maker(cs) = SigmaZetaMeshFromVMADCP(vmadcp);
    mesh_maker(cs).filter = ef(cs);
    mesh_maker(cs).bathymetry = bathy(cs);
    mesh_maker(cs).xs = xs(cs);
    mesh(cs) = mesh_maker(cs).get_mesh();
    vel_solver(cs) = TimeBasedVelocitySolver(mesh(cs),vmadcp,ef(cs),xs(cs),bathy(cs));
    vel{cs} = vel_solver(cs).get_velocity();
end

%%

hold on
for i = 1:numel(ef)
plot3(mesh(i),vel{i}(:,1));
bathy(i).plot;
colormap default
end
hold off

