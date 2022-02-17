%% Read, filter, plot ADCP data
% Created by Henk Jongbloed

clearvars % Clear all variables from the workspace
close all % Close all figures
clc
cd('C:\Users\jongb013\OneDrive - WageningenUR\PhDHenkJongbloed\2. Programming\Code\DataAnalysis')

addPaths
%% Import Data

% August 2014: Hartelkanaal (Verbeek)
load H_201408.mat %Imports waterlevels Spijkenisse: 1 Aug 2014 until 31 Aug 2014, 10 minute intervals
load ADCP_20140812.mat %1: HAK 2: OMS 3: OMN  , not needed (already imported in MSH_20140812)
load MSH_20140812.mat %1: HAK 2: OMS 3: OMN (only velocity data)

% September 2015: Nieuwe Waterweg (Poelman)
load H_2015091415.mat % Imports waterlevels Geulhaven, 14-Sep-15 00:00 until 15-Sep-15 12:24, 1 min interval
load CTD_20150914.mat  %See msh_CTD , not needed (already imported in MSH_20150914)
load MSH_20150914.mat %1: NM 2: OM 3: RWW (salt and velocity data)

% datstr = '20150914';

% msh = MSH_20140812;
% h = H_201408;

msh = MSH_20150914;
h = H_2015091415;

%% Pressure: Main Fig
% Create colormap
saltmap = brewermap(15, 'YlOrBr'); % Brewermap.m (available through Github, use the same code); map = brewermap(N,scheme)
% colormap(cmap)
clim(1) = min([min(msh{1}.CTD_int.Pressure_ex(:)), min(msh{2}.CTD_int.Pressure_ex(:)), min(msh{3}.CTD_int.Pressure_ex(:))]);
clim(2) = max([max(msh{1}.CTD_int.Pressure_ex(:)), max(msh{2}.CTD_int.Pressure_ex(:)), max(msh{3}.CTD_int.Pressure_ex(:))]);
saltFig = figure('units','normalized','outerposition',[0 0 1 1]);
for ct=1:3
    for hour=1:13
        sp = sub2ind([3,13],ct,hour);
        subplot(13,3,sp);
        if hour==1&&ct==1
            title('NM')
        elseif hour==1&&ct==2
            title('OM')
        elseif hour==1&&ct==3
            title('RWW')
        end
        patch(msh{ct}.p.N(:,:,1),...    % Plot salinity for N, Z
            msh{ct}.p.Z(:,:,1),...
            msh{ct}.CTD_int.Pressure_ex((msh{ct}.p.progfgood(:,hour)))', 'FaceColor', 'flat');
        caxis(clim);
        hold on
        plot((msh{ct}.p.nbed),...       % Plot channel bed
            msh{ct}.p.zbed,'-k', 'LineWidth', 1.15);
        hold off
        axis tight
        %         print(sp)
    end
end
sgtitle('Pressure: Measured')
colormap(saltmap);
cb = colorbar;

%% Pressure: Hydrostatic Approximation: Constant Density
% Create colormap
saltmap = brewermap(15, 'YlOrBr'); % Brewermap.m (available through Github, use the same code); map = brewermap(N,scheme)
% colormap(cmap)
clim(1) = min([min(msh{1}.CTD_int.Pressure_ex(:)), min(msh{2}.CTD_int.Pressure_ex(:)), min(msh{3}.CTD_int.Pressure_ex(:))]);
clim(2) = max([max(msh{1}.CTD_int.Pressure_ex(:)), max(msh{2}.CTD_int.Pressure_ex(:)), max(msh{3}.CTD_int.Pressure_ex(:))]);
saltFig = figure('units','normalized','outerposition',[0 0 1 1]);
for ct=1:3
    for hour=1:13
        sp = sub2ind([3,13],ct,hour);
        subplot(13,3,sp);
        if hour==1&&ct==1
            title('NM')
        elseif hour==1&&ct==2
            title('OM')
        elseif hour==1&&ct==3
            title('RWW')
        end
        patch(msh{ct}.p.N(:,:,1),...    % Plot salinity for N, Z
            msh{ct}.p.Z(:,:,1),...
            msh{ct}.p.Z(:,:,1), 'FaceColor', 'flat');
        caxis(clim);
        hold on
        plot((msh{ct}.p.nbed),...       % Plot channel bed
            msh{ct}.p.zbed,'-k', 'LineWidth', 1.15);
        hold off
        axis tight
        %         print(sp)
    end
end
sgtitle('Pressure: Measured')



colormap(saltmap);
cb = colorbar;