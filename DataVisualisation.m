%% Read, filter, plot ADCP data
% Created by Henk Jongbloed

clearvars % Clear all variables from the workspace
close all % Close all figures
clc

%% Import Data

% load data - water levels measured at station Spijkenisse (downloaded from
% https://waterinfo.rws.nl/#!/nav/index/)
% Time in UTC+1 (MET)
addPaths
load H_201408.mat   %Imports the struct h
load WURData_processed.mat      %Imports cell arrays msh and adcp
% load 20150914_h_geulhaven.mat
%% Unpack plotting data

% transect = 1; %1, 2 or 3
% In our data, this is the case:
CS = 3;% - number of cross-section
E = [length(adcp{1}.FileNumber), length(adcp{2}.FileNumber), length(adcp{3}.FileNumber)];% - number of ensembles
R = [size(msh{1}.pars,1), size(msh{2}.pars,1), size(msh{3}.pars,1)];%- maximum number of rows in output mesh
C = [size(msh{1}.pars,2), size(msh{2}.pars,2), size(msh{3}.pars,2)];%- number of columns in output mesh
T = [size(msh{1}.pars,3), size(msh{2}.pars,3), size(msh{3}.pars,3)];%- number of time steps in output
P = [size(msh{1}.pars,4), size(msh{2}.pars,4), size(msh{3}.pars,4)];%- number of fitted parameter (i.e. sum of the parameters for
                                %x,y and z components of the velocity)
Q = [size(msh{1}.p.X,2), size(msh{2}.p.X,2), size(msh{3}.p.X,2)];%- number of cells in the cross-section



%% Visualize data - Water Level

figure;
% subplot(3,1,1)
subplot(2,1,1)
hold on
plot(msh{1, 1}.eta,'k','LineWidth',2)
plot(msh{1, 1}.maxeta)
plot(msh{1, 1}.mineta)

% subplot(3,1,2)
% hold on
plot(msh{2, 1}.eta,'k--','LineWidth',2)
plot(msh{2, 1}.maxeta)
plot(msh{2, 1}.mineta)

% subplot(3,1,3)
% hold on
plot(msh{3, 1}.eta,'k-.','LineWidth',2)
plot(msh{3, 1}.maxeta)
plot(msh{3, 1}.mineta)
legend('Hartel Canal: Mean','Hartel Canal: Max','Hartel Canal: Min'...
    ,'Old Meuse South: Mean','Old Meuse South: Max','Old Meuse South: Min'...
    ,'Old Meuse North: Mean','Old Meuse North: Max','Old Meuse North: Min')

grid on
xlabel('time step [h]')
ylabel('water level [m]')
title('Tide in the three transects and at Spijkenisse')

subplot(2,1,2)
plot(h.level)  
xlabel('Days in 2014')
ylabel('water level [m]')
grid on


%% Visualize mesh and river bed
figure;
grid on
hold on
plot(msh{1}.p.nbed, msh{1}.p.zbed,'k')% add bed geometry %HartelCanal
plot(msh{2}.p.nbed, msh{2}.p.zbed,'k--')% add bed geometry %Old Meuse
plot(msh{3}.p.nbed, msh{3}.p.zbed,'k.-')% add bed geometry %Old Meuse
legend('Hartel Canal','Old Meuse South','Old Meuse North')
title('River beds')
xlabel('Cross-section')
ylabel('Depth')

%% Visualize data - Init

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

%% Example figure: Cell vertices
fig           = figure;
u            = fig.Units;
fig.Units    = 'centimeters';
% set(fig, 'Position', position);
% set(gcf, 'Position', get(0, 'Screensize'));
% cmap = brewermap(20, 'RdBu');
%%%%%%%%%%%%%%%%%%% Plots with Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ct = 1:3
    for t = 1
        fig          = subplot(3,1,ct);
        fig.Units    = 'centimeters';
%         set(fig, 'Position', [2+(ct-1)*(s+ax_h), position(4)-(1+t)*1*(s+ax_v)-3, ax_h, ax_v]);
        
        %%%%%%%%%%%%% Use function 'patch' to plot the grid with velocities
        p        = plot(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,t),'.');
        %%%%%%%%%%%%%
        hold on % This enable plotting something else 'on top of' the first plot.
        plot(msh{ct}.p.nbed, msh{ct}.p.zbed, '-k')% add bed geometry
%         plot(
%         caxis([-1,1]);              % Set the limits of velocity-values
%         set(p, 'edgecolor','none')  % Removes black lines of grid
        ylim([-15 2])
        xlim([-1 350])
%         colormap(cmap)
%         if ct==1 % Here, t=# is added for each row
%             handle=title(['\fontsize{8} t = ', num2str(t)],'FontWeight','Normal'); handle.Units='centimeters'; set(handle,'Position',[-1.5, 0.5, 1]);
            ylabel('Depth [m]')
%         else
%             set(gca,'YTickLabel',[]);
%         end
% %         
%         if t~=13
%             set(gca,'XTickLabel',[]);
%         else
            xlabel('x-coordinate [m]')
%         end
        %axis tight
        
%         if t==1 && ct==2
%         end
    end
end
sgtitle(['\fontsize{10} Cell vertices '],'FontWeight','Normal'); %handle.Units='centimeters'; set(handle,'Position',[1.3,1.6, 1]);

% Include legend:
% lgd       = colorbar;  
% lgd.Units = 'centimeters'; 
% set(lgd, 'Position', [2.2+(3)*(s+ax_h), position(4)-(3)*(s+ax_v)-3, 0.4, 2*ax_v]);
% lgd.Label.String = 'u [m s^{-1}]' ;
% 
% %export_fig([ROOTFOLDER, '\Figures\Velocity'], '-png', '-transparent', '-r464')
















%% Example figure: x-velocity
fig           = figure;
u            = fig.Units;
fig.Units    = 'centimeters';
set(fig, 'Position', position);
% set(gcf, 'Position', get(0, 'Screensize'));
cmap = brewermap(20, 'RdBu');
%%%%%%%%%%%%%%%%%%% Plots with Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ct = 1:3
    for t = 1:13
        fig          = subplot(no_figs1,no_figs2,36);
        fig.Units    = 'centimeters';
        set(fig, 'Position', [2+(ct-1)*(s+ax_h), position(4)-(1+t)*1*(s+ax_v)-3, ax_h, ax_v]);
        
        %%%%%%%%%%%%% Use function 'patch' to plot the grid with velocities
        p        = patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),...
            dir(ct)* msh{ct}.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),1))');
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
            handle=title(['\fontsize{10} Flow velocity in x-direction'],'FontWeight','Normal'); handle.Units='centimeters'; set(handle,'Position',[1.3,1.6, 1]);
        end
    end
end
% Include legend:
lgd       = colorbar;  
lgd.Units = 'centimeters'; 
set(lgd, 'Position', [2.2+(3)*(s+ax_h), position(4)-(3)*(s+ax_v)-3, 0.4, 2*ax_v]);
lgd.Label.String = 'u [m s^{-1}]' ;
% 
% %export_fig([ROOTFOLDER, '\Figures\Velocity'], '-png', '-transparent', '-r464')

%% Example figure: y-velocity
fig           = figure;
u            = fig.Units;
fig.Units    = 'centimeters';
set(fig, 'Position', position);
% set(gcf, 'Position', get(0, 'Screensize'));
cmap = brewermap(20, 'RdBu');
%%%%%%%%%%%%%%%%%%% Plots with Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ct = 1:3
    for t = 1:13
        fig          = subplot(no_figs1,no_figs2,36);
        fig.Units    = 'centimeters';
        set(fig, 'Position', [2+(ct-1)*(s+ax_h), position(4)-(1+t)*1*(s+ax_v)-3, ax_h, ax_v]);
        
        %%%%%%%%%%%%% Use function 'patch' to plot the grid with velocities
        p        = patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),...
            dir(ct)* msh{ct}.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),2))');
        %%%%%%%%%%%%
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
            handle=title(['\fontsize{10} Flow velocity in y-direction'],'FontWeight','Normal'); handle.Units='centimeters'; set(handle,'Position',[1.3,1.6, 1]);
        end
    end
end
% Include legend:
lgd       = colorbar;  
lgd.Units = 'centimeters'; 
set(lgd, 'Position', [2.2+(3)*(s+ax_h), position(4)-(3)*(s+ax_v)-3, 0.4, 2*ax_v]);
lgd.Label.String = 'u [m s^{-1}]' ;
% 
% %export_fig([ROOTFOLDER, '\Figures\Velocity'], '-png', '-transparent', '-r464')

%% Example figure: z-velocity
fig           = figure;
u            = fig.Units;
fig.Units    = 'centimeters';
set(fig, 'Position', position);
% set(gcf, 'Position', get(0, 'Screensize'));
cmap = brewermap(20, 'RdBu');
%%%%%%%%%%%%%%%%%%% Plots with Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ct = 1:3
    for t = 1:13
        fig          = subplot(no_figs1,no_figs2,36);
        fig.Units    = 'centimeters';
        set(fig, 'Position', [2+(ct-1)*(s+ax_h), position(4)-(1+t)*1*(s+ax_v)-3, ax_h, ax_v]);
        
        %%%%%%%%%%%%% Use function 'patch' to plot the grid with velocities
        p        = patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),...
            dir(ct)* msh{ct}.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),3))');
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
            handle=title(['\fontsize{10} Flow velocity in z-direction'],'FontWeight','Normal'); handle.Units='centimeters'; set(handle,'Position',[1.3,1.6, 1]);
        end
    end
end
% Include legend:
lgd       = colorbar;  
lgd.Units = 'centimeters'; 
set(lgd, 'Position', [2.2+(3)*(s+ax_h), position(4)-(3)*(s+ax_v)-3, 0.4, 2*ax_v]);
lgd.Label.String = 'u [m s^{-1}]' ;
% 
% %export_fig([ROOTFOLDER, '\Figures\Velocity'], '-png', '-transparent', '-r464')

%% Example figure: norm-velocity
fig           = figure;
u            = fig.Units;
fig.Units    = 'centimeters';
set(fig, 'Position', position);
% set(gcf, 'Position', get(0, 'Screensize'));
cmap = brewermap(20, 'RdBu');
%%%%%%%%%%%%%%%%%%% Plots with Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ct = 1:3
    for t = 1:13
        fig          = subplot(no_figs1,no_figs2,36);
        fig.Units    = 'centimeters';
        set(fig, 'Position', [2+(ct-1)*(s+ax_h), position(4)-(1+t)*1*(s+ax_v)-3, ax_h, ax_v]);
        
        %%%%%%%%%%%%% Use function 'patch' to plot the grid with velocities
        p        = patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),...
           sqrt((dir(ct)* msh{ct}.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),1))').^2+...
                (dir(ct)* msh{ct}.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),2))').^2+...
                (dir(ct)* msh{ct}.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),3))').^2));
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
            handle=title(['\fontsize{10} Flow velocity: norm'],'FontWeight','Normal'); handle.Units='centimeters'; set(handle,'Position',[1.3,1.6, 1]);
        end
    end
end
% Include legend:
lgd       = colorbar;  
lgd.Units = 'centimeters'; 
set(lgd, 'Position', [2.2+(3)*(s+ax_h), position(4)-(3)*(s+ax_v)-3, 0.4, 2*ax_v]);
lgd.Label.String = 'u [m s^{-1}]' ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Example figure: x-velocity N-dir
fig           = figure;
u            = fig.Units;
fig.Units    = 'centimeters';
set(fig, 'Position', position);
% set(gcf, 'Position', get(0, 'Screensize'));
cmap = brewermap(20, 'RdBu');
%%%%%%%%%%%%%%%%%%% Plots with Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ct = 1:3
    for t = 1:13
        fig          = subplot(no_figs1,no_figs2,36);
        fig.Units    = 'centimeters';
        set(fig, 'Position', [2+(ct-1)*(s+ax_h), position(4)-(1+t)*1*(s+ax_v)-3, ax_h, ax_v]);
        
        %%%%%%%%%%%%% Use function 'patch' to plot the grid with velocities
        p        = patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),...
            dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),1))');
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
            handle=title(['\fontsize{10} Flow velocity in nx-direction'],'FontWeight','Normal'); handle.Units='centimeters'; set(handle,'Position',[1.3,1.6, 1]);
        end
    end
end
% Include legend:
lgd       = colorbar;  
lgd.Units = 'centimeters'; 
set(lgd, 'Position', [2.2+(3)*(s+ax_h), position(4)-(3)*(s+ax_v)-3, 0.4, 2*ax_v]);
lgd.Label.String = 'u [m s^{-1}]' ;
% 
% %export_fig([ROOTFOLDER, '\Figures\Velocity'], '-png', '-transparent', '-r464')

%% Example figure: y-velocity N-dir
fig           = figure;
u            = fig.Units;
fig.Units    = 'centimeters';
set(fig, 'Position', position);
% set(gcf, 'Position', get(0, 'Screensize'));
cmap = brewermap(20, 'RdBu');
%%%%%%%%%%%%%%%%%%% Plots with Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ct = 1:3
    for t = 1:13
        fig          = subplot(no_figs1,no_figs2,36);
        fig.Units    = 'centimeters';
        set(fig, 'Position', [2+(ct-1)*(s+ax_h), position(4)-(1+t)*1*(s+ax_v)-3, ax_h, ax_v]);
        
        %%%%%%%%%%%%% Use function 'patch' to plot the grid with velocities
        p        = patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),...
            dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),2))');
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
            handle=title(['\fontsize{10} Flow velocity in ny-direction'],'FontWeight','Normal'); handle.Units='centimeters'; set(handle,'Position',[1.3,1.6, 1]);
        end
    end
end
% Include legend:
lgd       = colorbar;  
lgd.Units = 'centimeters'; 
set(lgd, 'Position', [2.2+(3)*(s+ax_h), position(4)-(3)*(s+ax_v)-3, 0.4, 2*ax_v]);
lgd.Label.String = 'u [m s^{-1}]' ;
% 
% %export_fig([ROOTFOLDER, '\Figures\Velocity'], '-png', '-transparent', '-r464')

%% Example figure: z-velocity N-dir
fig           = figure;
u            = fig.Units;
fig.Units    = 'centimeters';
set(fig, 'Position', position);
% set(gcf, 'Position', get(0, 'Screensize'));
cmap = brewermap(20, 'RdBu');
%%%%%%%%%%%%%%%%%%% Plots with Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ct = 1:3
    for t = 1:13
        fig          = subplot(no_figs1,no_figs2,36);
        fig.Units    = 'centimeters';
        set(fig, 'Position', [2+(ct-1)*(s+ax_h), position(4)-(1+t)*1*(s+ax_v)-3, ax_h, ax_v]);
        
        %%%%%%%%%%%%% Use function 'patch' to plot the grid with velocities
        p        = patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),...
            dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),3))');
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
            handle=title(['\fontsize{10} Flow velocity in nz-direction'],'FontWeight','Normal'); handle.Units='centimeters'; set(handle,'Position',[1.3,1.6, 1]);
        end
    end
end
% Include legend:
lgd       = colorbar;  
lgd.Units = 'centimeters'; 
set(lgd, 'Position', [2.2+(3)*(s+ax_h), position(4)-(3)*(s+ax_v)-3, 0.4, 2*ax_v]);
lgd.Label.String = 'u [m s^{-1}]' ;
% 
% %export_fig([ROOTFOLDER, '\Figures\Velocity'], '-png', '-transparent', '-r464')

%% Example figure: norm-velocity N-dir
fig           = figure;
u            = fig.Units;
fig.Units    = 'centimeters';
set(fig, 'Position', position);
% set(gcf, 'Position', get(0, 'Screensize'));
cmap = brewermap(20, 'RdBu');
%%%%%%%%%%%%%%%%%%% Plots with Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ct = 1:3
    for t = 1:13
        fig          = subplot(no_figs1,no_figs2,36);
        fig.Units    = 'centimeters';
        set(fig, 'Position', [2+(ct-1)*(s+ax_h), position(4)-(1+t)*1*(s+ax_v)-3, ax_h, ax_v]);
        
        %%%%%%%%%%%%% Use function 'patch' to plot the grid with velocities
        p        = patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),...
           sqrt((dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),1))').^2+...
                (dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),2))').^2+...
                (dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),3))').^2));
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
            handle=title(['\fontsize{10} Flow velocity: n-norm'],'FontWeight','Normal'); handle.Units='centimeters'; set(handle,'Position',[1.3,1.6, 1]);
        end
        
            diff(ct,t) = abs(norm(sqrt((dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),1))').^2+...
                (dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),2))').^2+...
                (dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),3))').^2))-norm(sqrt((dir(ct)* msh{ct}.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),1))').^2+...
                (dir(ct)* msh{ct}.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),2))').^2+...
                (dir(ct)* msh{ct}.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),3))').^2)));
        %Norm is preserved
    end
end
% Include legend:
lgd       = colorbar;  
lgd.Units = 'centimeters'; 
set(lgd, 'Position', [2.2+(3)*(s+ax_h), position(4)-(3)*(s+ax_v)-3, 0.4, 2*ax_v]);
lgd.Label.String = 'u [m s^{-1}]' ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


