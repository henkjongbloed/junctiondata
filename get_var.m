function [evvar, X, H] = get_var(var, res, lim)

X = get_coords(lim, res);
V = var.evaluate(1, T = X.T(:), N = X.Y(:), Sig = X.Sig(:), extrapolate=true);
evvar = cell([1,size(V,2)]);
for d=1:size(V,2)
    evvar{d} = reshape(V(:,d), res);
end
[Hf, ~, Zbf] = get_H(X, var, 1); % Own function to derive H(t,y) from X and adcptools variables.
H = repmat(Hf, [1,1,numel(X.sig)]);
X.Z = Zbf + X.Sig.*H;

end