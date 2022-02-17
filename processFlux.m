function F = processFlux(U, S)

[u, v, ~] = cells2mat(U, 'vel');


u =  decomposeVar(u', U);
v =  decomposeVar(v', U);

s = cells2mat(S, 'sal');

s = decomposeVar(s', U);

H = decomposeH(U);

[Lat, Loc] = SFD(H, u, v, s);

F.u = U; F.v = v; F.s = s; F.H = H; F.Lat = Lat; F.Loc = Loc;

end

