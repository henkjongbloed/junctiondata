function [Lat, Loc] = SFD(H, u, v, s)

Lat(:,1) = LateralDecomposition(H, u, s);
Lat(:,2) = LateralDecomposition(H, v, s);

Loc = Kjerfve(H, u, v, s);



% LatDec = B*LatDec;


% TA.FT = LatDec;
% TA.F = F;
% TA.u = xcell; TA.s = s; TA.v = v;


%[TA.A{1,1}, TA.A{1,2}] = getAmpPhase(xcell{2,2},t');
%[TA.A{2,1}, TA.A{2,2}] = getAmpPhase(v{2,2},t');
%[TA.A{3,1}, TA.A{3,2}] = getAmpPhase(s{2,2},t');
end