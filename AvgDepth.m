



% U = vel(:,1);
% V = vel(:,2);
% 
% 
% DA.uv = NaN(nr,nc); DA.vv = NaN(nr,nc); 
% DA.sv = NaN(nr,nc); DA.fv = NaN(nr,nc); 
% 
% 
% 
% for c = 1:nc
%     Uc = U(ci == c); %column c
%     Vc = V(ci == c); %column c
% 
%     DA.u0(c) = mean(Uc,'omitnan'); DA.uv(1:size(Uc),c) = Uc - DA.u0(c);
%     DA.v0(c) = mean(Vc,'omitnan'); DA.vv(1:size(Vc),c) = Vc - DA.v0(c);
% 
%     Sc = S(ci == c);
%     DA.s0(c) = mean(Sc,'omitnan'); DA.sv(1:size(Sc),c) = Sc - DA.s0(c);
%     Fc = F(ci == c);
%     DA.f0(c) = mean(Fc,'omitnan'); DA.fv(1:size(Fc),c) = Fc - DA.f0(c);
% end
% 
% DA.U = U; DA.S = S; DA.F = F;
% DA.V = V;
% DA.uCSA = mean(U,'omitnan'); DA.vCSA = mean(V,'omitnan'); 
% 
% DA.sCSA = mean(S,'omitnan'); DA.fCSA = mean(F,'omitnan');