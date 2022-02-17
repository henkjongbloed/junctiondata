function plotMSH(msh)

%% Velocity: Main Fig

dir           = [-1, 1, 1]; % Direction of flow
velmap = brewermap(20, 'RdBu');
velFig = figure('units','normalized','outerposition',[0 0 1 1]);
for ct=1:3
    for hour=1:13
        caxis([-1,1]);
        sp = sub2ind([3,13],ct,hour);
        subplot(13,3,sp);
        if hour==1&&ct==1
            title('NM')
        elseif hour==1&&ct==2
            title('OM')
        elseif hour==1&&ct==3
            title('RWW')
        end
        p        = patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),...
            dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood(:,msh{ct}.adcp_good(hour)))', 'FaceColor', 'flat','LineStyle','none');
        hold on
        plot((msh{ct}.p.nbed),...       % Plot channel bed
            msh{ct}.p.zbed,'-k', 'LineWidth', 1.15);
        hold off
        %         caxis([0,25]);
        axis tight
        %         print(sp)
    end
end
sgtitle('Velocity in o-direction')
colormap(velmap);

cb = colorbar; cb.Units = 'centimeters'; cb.Location = 'eastoutside';

%% Figure of salinity in cross-sections: Main Fig
% msh = MSH_20140812;

% Create colormap
saltmap = brewermap(15, 'YlOrBr'); % Brewermap.m (available through Github, use the same code); map = brewermap(N,scheme)
% colormap(cmap)
clim(1) = min([min(msh{1}.CTD_int.Salinity_ex(:)), min(msh{2}.CTD_int.Salinity_ex(:)), min(msh{3}.CTD_int.Salinity_ex(:))]);
clim(2) = max([max(msh{1}.CTD_int.Salinity_ex(:)), max(msh{2}.CTD_int.Salinity_ex(:)), max(msh{3}.CTD_int.Salinity_ex(:))]);

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
        %         patch(msh{ct}.p.N(:,:,1),...    % Plot salinity for N, Z
        %               msh{ct}.p.Z(:,:,1),...
        %               msh{ct}.CTD_int.Sal_ex((msh{ct}.p.progfgood(:,msh{ct}.adcp_good(hour))))', 'FaceColor', 'flat');
        patch(msh{ct}.p.N(:,:,1),...    % Plot salinity for N, Z
            msh{ct}.p.Z(:,:,hour),...
            msh{ct}.CTD_int.Salinity_ex((msh{ct}.p.progfgood(:,hour)))', 'FaceColor', 'flat');
        
        caxis([0,25]);
        hold on
        plot((msh{ct}.p.nbed),...       % Plot channel bed
            msh{ct}.p.zbed,'-k', 'LineWidth', 1.15);
        hold off
        axis tight
        %         print(sp)
    end
end
sgtitle('Salinity [psu]')
colormap(saltmap);
cb = colorbar;


%% Density: Main Fig
% Create colormap
saltmap = brewermap(15, 'YlOrBr'); % Brewermap.m (available through Github, use the same code); map = brewermap(N,scheme)
% colormap(cmap)
clim(1) = min([min(msh{1}.CTD_int.Density_ex(:)), min(msh{2}.CTD_int.Density_ex(:)), min(msh{3}.CTD_int.Density_ex(:))]);
clim(2) = max([max(msh{1}.CTD_int.Density_ex(:)), max(msh{2}.CTD_int.Density_ex(:)), max(msh{3}.CTD_int.Density_ex(:))]);
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
        %         patch(msh{ct}.p.N(:,:,1),...    % Plot salinity for N, Z
        %               msh{ct}.p.Z(:,:,1),...
        %               msh{ct}.CTD_int.Sal_ex((msh{ct}.p.progfgood(:,msh{ct}.adcp_good(hour))))', 'FaceColor', 'flat');
        patch(msh{ct}.p.N(:,:,1),...    % Plot salinity for N, Z
            msh{ct}.p.Z(:,:,1),...
            msh{ct}.CTD_int.Density_ex((msh{ct}.p.progfgood(:,hour)))', 'FaceColor', 'flat');
        
        caxis(clim);
        hold on
        plot((msh{ct}.p.nbed),...       % Plot channel bed
            msh{ct}.p.zbed,'-k', 'LineWidth', 1.15);
        hold off
        axis tight
        %         print(sp)
    end
end
sgtitle('Density [kg/m^3]')

colormap(saltmap);
cb = colorbar;

%% Temperature: Main Fig
% Create colormap
saltmap = brewermap(15, 'YlOrBr'); % Brewermap.m (available through Github, use the same code); map = brewermap(N,scheme)
% colormap(cmap)
clim(1) = min([min(msh{1}.CTD_int.Temperature_ex(:)), min(msh{2}.CTD_int.Temperature_ex(:)), min(msh{3}.CTD_int.Temperature_ex(:))]);
clim(2) = max([max(msh{1}.CTD_int.Temperature_ex(:)), max(msh{2}.CTD_int.Temperature_ex(:)), max(msh{3}.CTD_int.Temperature_ex(:))]);
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
        %         patch(msh{ct}.p.N(:,:,1),...    % Plot salinity for N, Z
        %               msh{ct}.p.Z(:,:,1),...
        %               msh{ct}.CTD_int.Sal_ex((msh{ct}.p.progfgood(:,msh{ct}.adcp_good(hour))))', 'FaceColor', 'flat');
        patch(msh{ct}.p.N(:,:,1),...    % Plot temperature for N, Z
            msh{ct}.p.Z(:,:,1),...
            msh{ct}.CTD_int.Temperature_ex((msh{ct}.p.progfgood(:,hour)))', 'FaceColor', 'flat');
        
        caxis(clim);
        hold on
        plot((msh{ct}.p.nbed),...       % Plot channel bed
            msh{ct}.p.zbed,'-k', 'LineWidth', 1.15);
        hold off
        axis tight
        %         print(sp)
    end
end
sgtitle('Temperature [C]')

colormap(saltmap);
cb = colorbar;

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
        %         patch(msh{ct}.p.N(:,:,1),...    % Plot salinity for N, Z
        %               msh{ct}.p.Z(:,:,1),...
        %               msh{ct}.CTD_int.Sal_ex((msh{ct}.p.progfgood(:,msh{ct}.adcp_good(hour))))', 'FaceColor', 'flat');
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
sgtitle('Pressure')

colormap(saltmap);
cb = colorbar;
