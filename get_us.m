function [u, s, D] = get_us(flow, salt, evres, reg_idx)
% evaluate u, s on a regular grid.



limu = {[0, flow.solver.model.periods(1,1)/(3600*24)];... %days (!)
    [min(flow.solver.mesh.n_left)+.5, max(flow.solver.mesh.n_right)-.5];...
    [0,1]};
% Flow

evresu = [evres(1), flow.solver.mesh.nverticals, flow.solver.mesh.max_ncells_vertical];
Xu = get_coords(limu, evresu);

u = get_var(flow, Xu, reg_idx);

% Salinity: Evaluate sigma from .1 up to .9
% limu
lims = limu;
lims{3} = [.1, .9];

evress = [evres(1), salt.solver.mesh.nverticals, salt.solver.mesh.max_ncells_vertical];
Xs19 = get_coords(lims, evress);

sc = get_var(salt, Xs19, reg_idx);

% Revert back to [0,1] sigma.
Xs = get_coords(limu, evress);
X = Xu;
% refine the salinity solution in the lateral using cubic interpolation.
% TODO - investigate result for multiple reg pars.
% s{1} = interpn(Xs.T,Xs.Y, Xs.Sig, sc{1}, X.T,X.Y, X.Sig);

s{1} = interpn(Xs.T, Xs.Y, Xs.Sig, sc{1}, X.T, X.Y, X.Sig, 'makima');

% Coordinate system of remaining analysis: X = Xu

[H, ~, Zb] = get_H(X, flow, 0);

if size(H, 3) == 1
    H = repmat(H, [1,1,numel(X.sig)]);
end
X.Z = Zb + X.Sig.*H;


D = Decomposition(X = X, H = H);

end