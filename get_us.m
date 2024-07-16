function [u, v, w, s, D] = get_us(flow, sal, evres)
M2T = flow.solver.model.periods(1,1)/(3600*24);

t = linspace(0, M2T, evres(1)); % in days
n = linspace(min(flow.solver.mesh.n_left)+.5, max(flow.solver.mesh.n_right)-.5, evres(2));
sig = linspace(0, 1, evres(3));

[T, N, Sig] = ndgrid(t,n,sig); Tq = T(:); Nq = N(:); Sigq = Sig(:);

U = flow.evaluate(1, T = Tq, N = Nq, Sig = Sigq, extrapolate=true);

u = reshape(U(:,1), evres);
v = reshape(U(:,2), evres);
w = reshape(U(:,3), evres);

X.t = t; X.y = n; X.sig = sig;
X.T = T; X.Y = N; X.Sig = Sig;

[Hf, ~, Zbf] = get_H(X, flow, 1); % Own function to derive H(t,y) from X and adcptools variables.

H = repmat(Hf, [1,1,numel(X.sig)]);
%Zb = repmat(Zbf, [1,1,numel(X.sig)]); 
X.Z = Zbf + X.Sig.*H;

S = sal.evaluate(1, T = Tq, N = Nq, Sig = Sigq, extrapolate=true);
S(S<0)=0;
s = reshape(S(:,1), evres);

D = Decomposition(X=X, H=H);
end