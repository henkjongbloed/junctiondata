%% Read, filter, plot ADCP data
% Created by Judith Poelman, 15-01-2019

clearvars % Clear all variables from the workspace
close all % Close all figures

% Add path to folder with ADCPTools
ROOTFOLDER = pwd;
% cd(ROOTFOLDER)
% startup_BScThesis(ROOTFOLDER)                                      % Calls function to add paths
addpath('./adcptools-code-r80-trunk')
addpath('./altmany-export_fig-dd9397c')
addpath('./Data')
addpath('./BrewerMap-master')

%% set preferences
nr_transects    = 3;     % [#] Number of transect in the dataset.
mesh_size_x     = 10;    % [m] Size of the grid: width of columns; can be changed
mesh_size_z     = 1.5;   % [m] Size of the grid: height of rows; can be changed

% load data - water levels measured at station Spijkenisse (downloaded from
% https://waterinfo.rws.nl/#!/nav/index/)
% Time in UTC+1 (MET)

% load 20140812_h_Spijkenisse.mat % h.date h.level

%% read in and filter raw adcp data (.000, .dat, .txt)
%
% adcp            = cell(3,1);                                    % Initialize empty cell array. In each cell the data from one transect will be stored.
%
% % ADCP data is read through function 'readDeployment'
% DEPNAME = 'OMHA'; % readDeployments read files starting with DEPNAME
%
% %load('adcp600.mat')
% for ct = 1:3
% PATH    = [ROOTFOLDER,'\Data\master',...
%     '\raai',num2str(ct), '\transect'];                          % reads files in PATH. PATH wordt gemaakt door een pad naar de folder waar de data staat
%                                                                 %(zou bij jou moeten eindigen in Data\Aquavision_ADCP_master_600.
%                                                                 %Vervolgens wordt '\raai + cijfer + \transect toegevoegd
% adcp{ct}        = readDeployment(DEPNAME,PATH); % Nu wordt alle data ingelezen, komt in workspace terecht.
% end

%% Interpolate data every hour to mesh 'msh'

% msh = cell(3,1); % Initialize empty cell array to store the results
%
% for ct = 1:nr_transects
%     time    = datenum(adcp{ct}.timeV);
%     eta     = interp1(h.date,h.level,time); % Calculate the waterlevel at each time step (using interpolation)
%     msh{ct}=procTrans2(...
%         adcp{ct},adcp{ct}.FileNumber,'ConventionalProcessing',false,...
%         'ConstantZetaMesh',true, 'TopMeshLowestEta', true,...
%         'DeltaN', mesh_size_x, 'DeltaZ', mesh_size_z, 'Eta',eta,...
%         'CumulateCrossings',false,'ShipReference','bt','EnableDebugging',false);
% end
%
% save('WURData_processed.mat', 'msh', 'adcp') % Can save the results for later use.

load 'WURData_processed.mat'

%% Example figure - initialize
tid_adcp = [1:13;    % Transect used to create a figure (for the second transect, there are some partial transect which are left out)
    1,2,3,6:15;
    1:13];
dir           = [-1, 1, 1]; % Direction of flow

position      = [1,0, 22, 30];        % Size and position of figure
no_figs1      = 6;                    % Total number of figures (no_figs1 times no_figs2)
no_figs2      = 6;
ax_v          = 1.5;                  % Size of subplot: height
ax_h          = 4;                    % Size of subplot: width
s             = 0.2;
%% Example figures
fig           = figure;
u            = fig.Units;
fig.Units    = 'centimeters';
set(fig, 'Position', position);

cmap = brewermap(20, 'RdBu');
%%%%%%%%%%%%%%%%%%% Plots with Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ct = 1:3
    for t = 1:13
        fig          = subplot(no_figs1,no_figs2,36);
        fig.Units    = 'centimeters';
        set(fig, 'Position', [2+(ct-1)*(s+ax_h), position(4)-(1+t)*1*(s+ax_v)-3, ax_h, ax_v]);
        
        %%%%%%%%%%%%% Use function 'patch' to plot the grid with velocities
        p        = patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood(:,tid_adcp(ct,t)))');
        %%%%%%%%%%%%%
        hold on % This enable plotting something else 'on top of' the first plot.
        plot(msh{ct}.p.nbed, msh{ct}.p.zbed, '-k')% add bed geometry
        caxis([-1,1]);              % Set the limits of velocity-values
        set(p, 'edgecolor','none')  % Removes black lines of grid
        ylim([-15 2])
        xlim([-1 350])
        colormap(cmap)
        if ct==1 % Here, t=# is added for each row
            handle=title(['\fontsize{8} t = ', num2str(t)],'FontWeight','Normal'); handle.Units='centimeters'; set(handle,'Position',[-1.5, 0.5, 1]);
            ylabel('Depth [m]')
        else
            set(gca,'YTickLabel',[]);
        end
        
        if t~=13
            set(gca,'XTickLabel',[]);
        else
            xlabel('x-coordinate [m]')
        end
        %axis tight
        
        if t==1 && ct==2
            handle=title(['\fontsize{10} Flow velocity'],'FontWeight','Normal'); handle.Units='centimeters'; set(handle,'Position',[1.3,1.6, 1]);
        end
    end
end
% Include legend:
lgd       = colorbar;  lgd.Units = 'centimeters'; set(lgd, 'Position', [2.2+(3)*(s+ax_h), position(4)-(3)*(s+ax_v)-3, 0.4, 2*ax_v]);
lgd.Label.String = 'u [m s^{-1}]' ;

%export_fig([ROOTFOLDER, '\Figures\Velocity'], '-png', '-transparent', '-r464')

