%% Read, filter, plot ADCP data
% Created by Henk Jongbloed

clearvars % Clear all variables from the workspace
close all % Close all figures
figsave = 0;
%clc
RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'WP2\TwoJunctions\code'])

% ds = 2; %Choose which dataset to analyze
% ds = input('Enter dataset'); % 1 XOR 2
% tak = input('Enter branch'); % 1 or 1:2 or 2:3 or 1:3 (all)

ds = 1;
tak = 1;

DS = {'OMHA14','NMOMNW15'}; DS = DS{ds};
BN = {{' Hartel Canal'; ' Old Meuse South'; ' Old Meuse North'},...
    {' New Meuse'; ' Old Meuse'; ' Rotterdam Waterway'}}; BN = BN{ds};

%% Import preprocessed semi-raw data (ADCP / CTD / H structs)

[adcp, ctd, figs] = import_data(RF, DS);

%% Data analysis
for j = 1:length(tak)
    b = tak(j);
    [V{j}, U{j}] = processVel(adcp{b}, figs, BN{b}); % V = VMADCP objects, U = generic velocity data (plus more)
    %     S{j} = processSal(U{j}, ctd{b}, DS); %S = Salinity data
    %     F{j} = processHours(U{j}, S{j});
    %     F{j} = processFlux(F{j});
    



%plot(U{1, 1}.T.velocity_model.names, U{1,1}.I)


% [M1, b1, Mj1, bj1] = GD2M(gather_dat);

% condest(M1)
% condest(M1'*M1)

% for j = 1:length(Mj)
%     frm(j) = norm(Mj{j} - Mj1{j}, 'fro');
%     frb(j) = norm(bj{j} - bj1{j}, 'fro');
% end
% plot_geom_all(U)
%% Plotting
plt=1;
if plt
    tak=1;
    %     plotCovParsx(U{i})
    plotPars_all(U{j}, {'New'})
end
%%
% figure;
% plot(cve)
% figure;
% plot(p)

%
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
% %
% % %% Tidal Ellipses
% % %
% plotTidalEllipses(U, F)
% plotTidalArrows(U, F)
%
% % plotTZ(U,F, -.6)
% plotTZ(U,F, 0)
% % plotTZ(U,F, .6)
%
% Save all figs automatically
figsave = 0;
if figsave
    figs =  findobj('type','figure');
    nfig = length(figs);
    FF = 'C:\Users\jongb013\Documents\PHD\2-Programming\WP2\TwoJunctions\figs\21-10\'; %RootFolder
    %cd([RF,'WP2\TwoJunctions\code'])
    % cd FF
    for j = 1:nfig
        f = gcf;
        st = f.Children.Title.String;
        dest = fullfile(FF,[num2str(j), '.png']);
        %         set(gcf, 'Position', get(0, 'Screensize'));
        saveas(f, dest);
        close;
    end
end
%
%
