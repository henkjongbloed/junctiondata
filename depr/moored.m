% moored.m

RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'WP2\TwoJunctions\code'])


addpath(genpath(strcat(RF,'Tools\adcptoolsGit'))); %path to ADCPTools
addpath(strcat(RF,'Tools\loess-master')); %path to loess-master
% addpath(strcat(RF,'Tools\T_Tide'))
addpath(strcat(RF,'Tools\BrewerMap-master')); %path to BrewerMap (see github)
% addpath(strcat(RF,'Tools\subaxis'));

F2 = 'WP2\TwoFrames\';
addpath(strcat(RF,'Tools\plotting\plotting'));

    % load new data
datloc = strcat(RF,F2, 'WUR\');
cd(datloc)
adcp{1} = rdi.readADCP('RWS1_000.000');
adcp{2} = rdi.readADCP('WUR1_000.000');