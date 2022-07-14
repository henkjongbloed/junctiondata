%% Read, filter, plot ADCP data
% Created by Henk Jongbloed

clearvars % Clear all variables from the workspace
close all % Close all figures
figsave = 0;
%clc
RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'WP2\TwoJunctions\code'])

% ds = 2; %Choose which dataset to analyze
ds = input('Enter dataset');
tak = input('Enter branch');
DS = {'OMHA14','NMOMNW15'}; DS = DS{ds}; 
BN = {{' Hartel Canal'; ' Old Meuse South'; ' Old Meuse North'},...
    {' New Meuse'; ' Old Meuse'; ' Rotterdam Waterway'}}; BN = BN{ds};

%% Import preprocessed semi-raw data (ADCP / CTD / H structs)

[adcp, ctd, figs] = import_data(RF, DS);

%% Data analysis
for j = 1:length(tak)
    b = tak(j);
    [V{j}, U{j}] = processVel(adcp{b}, figs, BN{b}); % V = VMADCP objects, U = generic velocity data (plus more)
    S{j} = processSal(U{j}, ctd{b}, DS); %S = Salinity data
    %F{j} = processHours(U{j}, S{j});
    %F{j} = processFlux(F{j});
end

%% Plotting
for i=1:length(tak)
    plotParsx(U{i})
end


% plotMap(U,S);
% 
% plotDA_Flux(tak, U, F);
% plotSubtidalFlow(tak, U);
% plotSubtidalSalt(tak, F, U);
% plot_geom_all(U, S)
% % plotPars(tak, U);
% plotT = 1:2:13;
% plotTimeUS(tak, plotT, F, U);
% plotCSA_Flux(tak, F, U);
% plotWL_CSA(F, tak, BN);
% 
% %% Tidal Ellipses
% % 
% plotTidalEllipses(U, F)
% plotTidalArrows(U, F)

% plotTZ(U,F, -.6)
%plotTZ(U,F, 0)
% plotTZ(U,F, .6)

%% Save all figs automatically
figsave = 0;
if figsave
    figs =  findobj('type','figure');
    nfig = length(figs);
    FF = 'C:\Users\jongb013\Documents\PHD\2-Programming\WP2\TwoJunctions\figs\26-4\'; %RootFolder
    %cd([RF,'WP2\TwoJunctions\code'])
    % cd FF
    for i = 1:nfig
        dest = fullfile(FF,[DS, num2str(i), '.png']);
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf, dest);
        close;
    end
end


