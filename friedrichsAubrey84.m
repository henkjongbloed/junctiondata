%FriedrichsAubrey1984


a = [1 0]; th = [pi pi];
v = [1 0]; ph = [0 0];

t = 1:1/60:24;

T2 = 12.42; O2 = 2*pi/T2;

A1 = a(1).*cos(O2.*t - th(1)) ;
A2 = a(2).*cos(2*O2.*t - th(2));
V1 = v(1).*cos(O2.*t - ph(1)) ;
V2 = v(2).*cos(2*O2.*t - ph(2));

figure;
subplot(2,1,1)
hold on
plot(t, A1, 'k:')
plot(t, A2, 'k--')
plot(t, A1+A2,'b-', 'LineWidth', 2)

title(['Surface: ', 'Ratios ', num2str(2*th(1) - th(2)), ' and ', num2str(a(2)/a(1))])

subplot(2,1,2)
hold on
plot(t, V1,'k:')
plot(t, V2,'k--')
plot(t, V1+V2,'b-', 'LineWidth', 2)
title(['Velocity: ', 'Ratios ', num2str(2*ph(1) - ph(2)), ' and ', num2str(v(2)/v(1))])

