% Analyze water levels

clearvars % Clear all variables from the workspace
close all % Close all figures
%clc
RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'WP2\TwoJunctions\code'])

DS = {'Maasvlakte', 'Suurhoffbrug', 'Harmsenbrug', 'Spijkenisse', 'Goidschalxoord', 'HoekvanHolland', 'Maassluis', 'Geulhaven', 'Vlaardingen', 'Rotterdam'};
DS = {'Spijkenisse', 'Geulhaven'};

H1 = import_waterlevels(RF, DS);

% dur = [datenum('11-Aug-2014 00:00:00'), datenum('13-Aug-2014 23:50:00')];
dur = [datenum('14-Sep-2015 00:00:00'), datenum('14-Sep-2015 00:00:00')];

%% Plotting
figure;

for i = 1:length(DS)
    subplot(length(DS),3,3*i-2)
    plot(H1{i}.t, H1{i}.wl)
    hold on
    plot([dur(1), dur(1)], [min(H1{i}.wl), max(H1{i}.wl)], 'y')
    plot([dur(2), dur(2)], [min(H1{i}.wl), max(H1{i}.wl)], 'y')
    xlabel('t')
    ylabel('m NAP')
    title(H1{i}.name)
    axis tight

    subplot(length(DS),3,3*i-1)
    idx = (H1{i}.t>dur(1) & H1{i}.t<dur(2));
    plot(H1{i}.t(idx), H1{i}.wl(idx))
    xlabel('t')
    ylabel('m NAP')
    title(H1{i}.name)
    axis tight
    
    Ts = 60*10; Tnyq = 2*Ts; Tpass = [2*365*24*3600, .1*24*3600];
    fs = 1./Ts;
    fNyquist = 1./fs;
    fpass = 1./Tpass;
    subplot(length(DS),3,3*i)
%     idx = (H1{i}.t>dur(1) & H1{i}.t<dur(2));
    fwl{i} = bandpass(H1{i}.wl, fpass, fs);
    plot(H1{i}.t, fwl{i})
    xlabel('t')
    ylabel('m NAP')
    title(H1{i}.name)
    axis tight
end

%% Data SP
% [nameu,fu,tidecon,xout] = t_tide(H2{6}.wl,1/6,H1{4}.t(1));
% t_tide(H1{i}.wl)

