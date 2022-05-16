function F = processHours(U, S)
% This function inter / extrapolates salinity values on a mesh that varies
% in time, such that we get (u, s) on all repeat transects, all mesh cells

hours = unique(U.adcpH);


mean_time = datenum(U.mesh_mean.time);
dt = hours - mean_time;
eta = U.eta(hours);

mesh = U.mesh_mean.mesh_at_water_level(eta); % compute expanded meshes

close all
for rt = 1:length(mesh)
%     mesh(rt).time = datetime(datestr(hours(rt)));
    SI = scatteredInterpolant(S.n(rt == (S.i+1)), S.sig(rt == (S.i+1)), S.S(rt == (S.i+1)), 'linear', 'nearest'); % Interpolator n, sigma, salt
    mn = mesh(rt).n_middle(mesh(rt).col_to_cell)'; 
    msig = mesh(rt).sig_center;
%     cell_centers = get_cell_centers();% center of cells in n, sigma coordinates
    mS = SI(mn, msig);
%     hold on
%     plot(mS)
    F.s{rt} = mS;
    F.u{rt} = U.T.velocity_model.get_velocity(U.pars{1, 1}, U.cov_pars{1,1}, 3600*24*dt(rt));
%     [vel, cov_vel] = T.velocity_model.get_velocity(U.pars{1,1}, U.cov_pars{1,1});
end

F.mesh = mesh;
F.mesh_mean = U.mesh_mean;
F.eta = U.eta(hours);
F.time = hours;
% F.dt = dt;

end




