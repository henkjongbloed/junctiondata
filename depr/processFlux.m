function F = processFlux(F)

[u, v, w] = cells2mat(F.u);
s = cells2mat(F.s);

u =  decomposeVar(u', F);
v =  decomposeVar(v', F);
w =  decomposeVar(w', F);



s = decomposeVar(s', F);

H = decomposeH(F);

% [Lat, Loc] = SFD(H, u, v, s);

Loc = Kjerfve(H, u, v, s);
Lat(:,1) = LateralDecomposition(H, u, s);
Lat(:,2) = LateralDecomposition(H, v, s);

F.uc = u; F.vc = v; F.wc = w;

F.sc = s; F.Hc = H; 

F.Lat = Lat; F.Loc = Loc;


XTL = {'F_{00}', 'F_{01}^*', 'F_{02}^*', 'F_{03}^*', 'F_{04}^*',...
    'F_{10}', 'F_{11}^*',...
    'F_{20}', 'F_{21}^*',...
    'F_{30}', 'F_{31}^*', 'F_{32}^*', 'F_{33}^*', 'F_{34}^*', 'F_{35}', 'F_{36}^*',...
    'F_{40}', 'F_{41}', 'F_{42}'};

Lat = Lat(:,1);

LatT = zeros(9, 1);
LatT([1, 3, 5, 7]) = Lat([1, 6, 8, 10]);
LatT(2) = sum(Lat(2:5));
LatT(4) = Lat(7);
LatT(6) = Lat(9);
LatT(7) = LatT(7) + Lat(15);
LatT(8) = sum(Lat([11:14, 16]));
LatT(9) = sum(Lat(17:19));

F.LatT = LatT;


LatAbs = zeros(9, 1);
LatAbs([1, 3, 5, 7]) = abs(Lat([1, 6, 8, 10]));
LatAbs(2) = sum(abs(Lat(2:5)));
LatAbs(4) = abs(Lat(7));
LatAbs(6) = abs(Lat(9));
LatAbs(7) = abs(LatAbs(7)) + abs(Lat(15));
LatAbs(8) = sum(abs(Lat([11:14, 16])));
LatAbs(9) = sum(abs(Lat(17:19)));

F.LatAbs = LatAbs;

F.LatDuoAbs(1) = sum(LatAbs([1, 3, 5, 7, 9])) ;
F.LatDuoAbs(2) = sum(LatAbs([2, 4, 6, 8])) ;

F.LatDuo(1) = sum(LatT([1, 3, 5, 7, 9])) ;
F.LatDuo(2) = sum(LatT([2, 4, 6, 8])) ;

end

