% coordinateTransform
clearvars; close all;
B = 300;
% y = 0:(2*B);

zb = @(y) 10./B.^2.*(y-B).^2; %max zb is 10, B must be in workspace
%trapz(y,zb(y))
eta = @(t) (15 + 2.*sin(2*pi./13.*t));

h = @(y,t) eta(t) - zb(y);
A = @(t) (2.*(eta(t)) - 20/3).*B;

sigma = @(z,y,t)  (z - zb(y))./(eta(t) - zb(y));

lambda = @(y,t) eta(t).*y./A(t) - 10./A(t)./3./B.^2.*((y-B).^3+B.^3);



% Continue here!
% y = 
% z = @(s,l,t) s.*eta(t) - zb(y
% 

%lambda = @(y) = 
%plot(y,zb)
%hold on
%plot(y,eta)