function [u, s, D] = get_us(flow, sal, evres)

lim = {[0, flow.solver.model.periods(1,1)/(3600*24)];...
    [min(flow.solver.mesh.n_left)+.5, max(flow.solver.mesh.n_right)-.5];...
    [0,1]};
% Flow
[u,X,H] = get_var(flow, evres, lim);

% Salinity: Evaluate sigma from .1 up to .9
lim{3} = [.0, .9];
[s,~, ~] = get_var(sal, evres, lim);


D = Decomposition(X=X, H=H);
end