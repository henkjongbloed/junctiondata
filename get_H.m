function [H, Wl, Zvec] = get_H(X, mesh, bathy, water_level, xs_proj)

% Gets the local, time-varying water depth from X and adcptools variables
% Wrapper for Solution.get_H


if xs_proj
    [xvec, yvec] = mesh.xs.sn2xy(zeros(size(X.y)), X.y); % convert to x,y -> TODO Enter correct s,n coord!!
else
    xvec = X.x;
    yvec = X.y; % convert to x,y -> TODO Enter correct s,n coord!!
end



zvec = bathy.get_bed_elev(xvec, yvec);  % get (time-indep) bathy
wl = water_level.get_water_level_model(datetime(X.t, 'ConvertFrom', 'datenum'));


[Wl, Zvec] = ndgrid(wl, zvec);

H = Wl - Zvec;

end