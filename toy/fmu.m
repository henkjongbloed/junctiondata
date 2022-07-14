function [u, R] = fmu(p, dr)

% creates a simple realization of u = [u,v,w]' at point r = r0 + dr
lp = length(p)/3;
lr = size(dr,2);
if lp == 1  %  only mean flow
    r = ones([lr,1]);
elseif lp == 4 %  first order approximation
    r = [ones([lr,1]), dr'];
elseif lp == 10  %  second order approximation
    r = [ones([lr,1]), dr', dr(1,:)'.^2/2, dr(1,:)'.*dr(2,:)', dr(1,:)'.*dr(3,:)', dr(2,:)'.^2/2, dr(2,:)'.*dr(3,:)', dr(3,:)'.^2];
end

R = zeros([3*lr, 3*lp]);
R(1:3:end, 1:lp) = r;
R(2:3:end, (lp+1):2*lp) = r;
R(3:3:end, (2*lp+1):3*lp) = r;

u = R*p;
% u = zeros([3,size(r,2)]);
% 
% 
% u(1,:) = p(1) + p(4)*r(2,:); % u varies linearly with x
% u(2,:) = p(2) + p(5)*r(3,:); % v varies linearly with y
% u(3,:) = p(3) + p(6)*r(4,:); % w varies linearly with z 

end