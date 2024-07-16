function [H, Wl, Zvec] = get_H(X, sol, xs_proj)


if xs_proj
    [xvec, yvec] = sol.solver.mesh.xs.sn2xy(zeros(size(X.y)), X.y); % convert to x,y -> TODO Enter correct s,n coord!!
else
    xvec = X.x;
    yvec = X.y; % convert to x,y -> TODO Enter correct s,n coord!!
end

zvec = sol.solver.bathy.get_bed_elev(xvec, yvec);  % get (time-indep) bathy
wl = sol.solver.bathy.water_level.get_water_level_model(datetime(X.t, 'ConvertFrom', 'datenum'));

[Wl, Zvec] = ndgrid(wl, zvec);

H = Wl - Zvec;

end