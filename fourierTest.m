%FOURIER TEST
clearvars

Fs = 100;           % Sampling frequency
t = 1:1/Fs:13;  % Time vector 
L = length(t);      % Signal length


T2 = 12.42; O2 = 2*pi/T2;
T4 = T2/2;  O4 = 2*pi/T4;

A2 = 1.5;
A4 = 1;

P2 = .5;
P4 = 1.5;

%t = 1:13;

u = A2*cos(O2*t - P2) + A4*cos(O4*t - P4);
u2 = u + .3*randn(1,length(u));

plot(t,u)
hold on
plot(t,u2)
% uf = fft(u);
% 
% f = Fs*(0:(L/2))/L;
% Tvec = 2*pi./f;
% 
% P = abs(abs(uf/L).^2);
% 
% 
% subplot(211)
% plot(t,u)
% subplot(212)
% plot(f,P(1:L/2+1)) 
%plot()

%This was the FFT approach, we now employ the normal equations.

A = [cos(O2*t)' sin(O2*t)' cos(O4*t)' sin(O4*t)'];
cond(A)
p = (A'*A)\A'*[u' u2'];

a(1,:) = sqrt(p(1,:).^2+p(2,:).^2);
a(2,:) = sqrt(p(3,:).^2+p(4,:).^2);

phi(1,:) = atan(p(2,:)./p(1,:));
phi(2,:) = atan(p(4,:)./p(3,:));

plot(t, a(1,1)*cos(O2*t - phi(1,1)) + a(2,1)*cos(O4*t - phi(2,1)));
plot(t, a(1,2)*cos(O2*t - phi(1,2)) + a(2,2)*cos(O4*t - phi(2,2)));

legend('Original', 'Perturbed', 'Est. Original', 'Est. Perturbed')


