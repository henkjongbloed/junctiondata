%modelmain
clc;
close all;
clearvars;

%% Generate artificial data

% r = [t, x, y, z] Coordinates of measured velocities, t in hours, rest in
% meters
n = 3000;
t = linspace(0, 13.5, n); %hours
r = [2*rand([1,n])-1; sin(2*pi*t)/2+1/2*(rand([1,n])-1); rand([1,n])-1]; %x, y, z (normalized): Measurement locations.
% figure
% scatter3(r(1,:), r(2,:), r(3,:))
r0 = mean(r,2);
dr = r - r0;
%% Forward model

%input parameter vector: 30 parameters
% A00 = [1,.1, -.2]';                   % u, v, w mean (size of vector = 3)
% A01 = rand([9,1]);    %[1, 0, 0, 0, -1, 0, 1, 0, 1]';  % u, v, w linearly varying (size of vector = 9)
% A02 = .1*rand([18,1]);               % u, v, w quadratically varying (size of vector = 18)
Ntide = 1;
Nspace = 1;

p = assemblep(Ntide, Nspace);
M = assembleM(t, dr, Ntide, Nspace);

%% ADCP measures radial velocities.
rb = [.2*(2*rand([1,n])-1); sin(2*pi*t)/2+1/6*(rand([1,n])-1/2); 0*rand([1,n])]; %x, y, z (normalized): Boat locations.

Mb = assembleMb(M, r, rb);

b = Mb*p; %Forward model: Radial velocity measurements
% [u, R] = fmu(p, dr);
% cond(Mb)

%% Inverse model direct

for i = 1:10
    b_noise = b + .05*(i-1)*randn(size(b));
    pest(:,i) = Mb\b_noise;
    er(i) = norm(pest(:,i)-p)./norm(p);
end
% s = fms(r);
figure;
subplot(4,1,1:3)
plot(pest)
hold on

plot(p,'r*', 'LineWidth', 3)
legend(strsplit(num2str(0:0.05:0.45)))
title('Parameter estimation unregularized')
subplot(4,1,4)

plot(0:0.05:0.45, er)
title('Relative error')

% To do: Temporal aspect: Tides
%% Inverse model Tikhonov
figure;
llist = [0, .5, 1, 2, 4, 10];
for j = 1:length(llist)
    for i = 1:10
        b_noise = b + .05*(i-1)*randn(size(b));
        lambda = llist(j);
        L = lambda*eye(size(Mb,2));
        pest(:,i)= (Mb'*Mb + L'*L)\Mb'*b_noise;
        er(i) = norm(pest(:,i)-p)./norm(p);
    end
    subplot(length(llist), 2, 2*j-1)
    plot(pest)
    hold on
    plot(p,'r*', 'LineWidth', 3)
    %     legend('pest','p')
    title(strcat('Parameter estimation Tikhonov: \lambda = ', num2str(llist(j))))
    subplot(length(llist), 2, 2*j)
    plot(0:0.05:0.45, er)
    title('Relative error')

end
% s = fms(r);

% To do: Temporal aspect: Tides

%% Inverse model SVD
[U, S, V] = svd(Mb, 'econ');
llist = round(size(Mb,2)*[.1, .2, .4, .6, .8, 1]);
figure;
for j = 1:length(llist)
    for i = 1:10
        b_noise = b + .05*(i-1)*randn(size(b));
        r = llist(j);%size(S,2);
        %     L = lambda*eye(size(Mb,2));
        pest(:,i)= V(:,1:r)*inv(S(1:r, 1:r))*U(:,1:r)'*b_noise;
        er(i) = norm(pest(:,i)-p)./norm(p);
    end
    subplot(length(llist), 2, 2*j-1)
    plot(pest)
    hold on
    plot(p,'r*', 'LineWidth', 3)
    %     legend('pest','p')
    title(strcat('Parameter estimation SVD: r = ', num2str(llist(j))))
    subplot(length(llist), 2, 2*j)
    plot(0:0.05:0.45, er)
    title('Relative error')

end
% s = fms(r);

% To do: Temporal aspect: Tides