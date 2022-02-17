%test AmpPhase
clearvars
%Time

T2 = 12.42; O2 = 2*pi/T2;
T4 = T2/2;  O4 = 2*pi/T4;
B = 300;
t = linspace(0, T2, 13)';
y = linspace(-B, B, 50);

A2 = @(y) 1.5 + y./(2.*B); P2 = @(y) pi/6 ;
A4 = @(y) .5; P4 = @(y) pi/4 + y./(4.*B);

u = @(t,y) A2(y).*cos(O2.*t + P2(y)) + A4(y).*cos(O4.*t + P4(y));

[T,Y] = meshgrid(t,y);

U = u(T,Y);

s = .3;

U = U + s.*randn(size(U));
%surf(T,Y,U)

[a, p] = getAmpPhase(U', t);

figure;
subplot(221)
hold on
plot(y, A2(y).*ones(size(y)));
plot(y, a(1,:));

subplot(222)
hold on

plot(y, A4(y).*ones(size(y)));
plot(y, a(2,:));

subplot(223)
hold on

plot(y, P2(y).*ones(size(y)));
plot(y, p(1,:));

subplot(224)
hold on

plot(y, P4(y).*ones(size(y)));
plot(y, p(2,:));

%u2 = @(t,y) A2(y).*cos(O2.*t + P2(y)) + A4(y).*cos(O4.*t + P4(y));



