%Tidal Analysis

% load hyperbolicChirp
% figure
% plot(t,hyperbolchirp)
% axis tight, xlabel('Seconds')
% ylabel('Amplitude')
% title('Hyperbolic Chirp')
% 
% Fs = 1/mean(diff(t));
% figure
% [cfs,f] = cwt(hyperbolchirp,Fs);
% cwt(hyperbolchirp,Fs);
% helperHyperbolicChirpPlot(cfs,f,t)
close all
clearvars
clc

addpath('./processedData')
file = 'test';
if strcmp(file,'H_201408')
    load(file)
    eta = H_201408.level;
    time = H_201408.date-H_201408.date(1); %one day = 1. (1 minute increments)
    loc = H_201408.location;
    t1 = H_201408.datestr{1, 1}  ;
    t2 = H_201408.datestr{end, 1}  ;
elseif strcmp(file,'H_2015091415')
    load(file)
    eta = H_2015091415.level;
    time = H_2015091415.date-H_2015091415.date(1); %one day = 1, so one hour = 1/24. (10 minute increments)
    loc = H_2015091415.location;
    t1 = H_2015091415.datestr{1, 1}  ;
    t2 = H_2015091415.datestr{end, 1}  ;
else %synthetic data
    days = 365;
    n = 144*365;
    time = linspace(0,365,n)';
    TM2 = (12+25.2/60)/24;
    TSN = 15;
    TY = 2*365/4;
    eta = 2*sin(2*pi/TM2*time);%+0.5*sin(2*pi/TSN*time) + 2*sin(2*pi/TY*time);
    t1 = 't=0'; t2 = 't=365'; loc = 'Timboektoe';
end
fs = round(1/(time(2)-time(1))); %times per day
tft = fft(eta);

n = length(eta);
power = abs(tft).^2/n;
%%
tft0 = fftshift(tft);         % shift y values
f0 = (-n/2:n/2-1)*(fs/n); % 0-centered frequency range
power0 = abs(tft0).^2/n;    % 0-centered power
ind = 1:length(f0);%find(abs(f0)<5);


%%
wl = 'bump';
fb = cwtfilterbank('SignalLength',n,'Wavelet',wl,'SamplingFrequency',fs);%,'FrequencyLimits',[0 6]); %freqz(fb)
% fbs = struct(fb);
[twt,tfreq,coi] = wt(fb,eta);
% [twt,tfreq,coi] = cwt(eta);

% figure
% pcolor(time,tfreq,abs(twt))
% shading flat
% % set(gca,'YScale','log')
% hold on
% plot(time,coi,'w-','LineWidth',1)
% axis tight
% xlabel('Time [day]')
% ylabel('Frequency [day^{-1}]')
% title(strcat('Wavelet Transform of tides at', 'Spijkenisse (August 2014)'))

%% Inverse
etawt = icwt(twt,wl,'SignalMean',mean(eta));
etaft = ifft(tft);

%% Plotting
figure
sgtitle(sprintf('Analysis of tide-influenced water levels in %s from %s until %s',loc,t1,t2))
ploty = 4;
plotx = 2;

subplot(4,2,1:2)
plot(time,eta)
grid on
axis tight
title('Water level data')
subplot(4,2,3)
plot(f0(ind),power0(ind))
grid on
axis tight
title('Fourier Power Spectrum')
subplot(4,2,5)
plot(time,etaft)
grid on
axis tight
title('Fourier Reconstruction')
subplot(4,2,7)
plot(time,etaft-eta)
grid on
axis tight
title('Fourier Reconstruction error')
subplot(4,2,4)
pcolor(time,tfreq,abs(twt))
shading flat
% % set(gca,'YScale','log')
hold on
plot(time,coi,'w-','LineWidth',1)
axis tight
hold off
title('Wavelet Transform')

subplot(4,2,6)
plot(time,etawt')
grid on
axis tight
title('Wavelet Reconstruction')
subplot(4,2,8)
plot(time,etawt'-eta)
grid on
axis tight
title('Wavelet Reconstruction Error')

