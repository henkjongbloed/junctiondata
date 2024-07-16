clearvars
close all
clc
addPaths
load MSH_20150914.mat 
load h.mat
load ctd.mat
                                                    % only time values

%% Interpolate salt profiles per hour
%figure;
msh = MSH_20150914;
% Start transect loop
for ct =1:3
    % start hour (time) loop
    for hour = 1:13
%         X=[];Y=[];T=[];Z=[];S=[];    % Initiate variables
        
        % for every hour and transect find indices
        ind        = find(and(ctd.Hour(:) == hour, ctd.Trans(:) == ct)); %Aangepast
        X          =  ctd.X(ind);
        Y          =  ctd.Y(ind);
        Z          =  ctd.Z_nap(ind);
        S          =  ctd.S(ind);
       
        %%%%%%%%% Step 1.  Convert transect to to n,z-coordinates %%%%%%%%%
        % n - along the transect, z - upwards
        [x,y]      = wgs2utm(Y,X);
        
        P_tmp      = [x,y,Z];                            % create CTD position data vector
         
        % Transformation of the coordinates to the mesh; following the method that was applied to the ADCP data (in procTrans2.m by B. Vermeulen) 
        P          = msh{ct}.Tvec(1)*(P_tmp(:,1)-msh{ct}.Pm(1))+msh{ct}.Tvec(2)*(P_tmp(:,2)-msh{ct}.Pm(2));
        P          = [P, P_tmp(:,3)];

        %%%%%%%%% Step 2. Interpolate data using scattered interpolant %%%%
        % Define mesh location at which interpolant will be evaluated

        nq         = msh{ct}.N(:,:,1)-max(max(msh{ct}.N(:,:,1)))./2; % N values in msh, transformed as Pm should be centre of the mesh. 
        zq         = msh{ct}.Z(:,:,1);
        
        % Define and plot interpolant ------------------------------------------
        f          = 1/50;%1/30;% Why 30?

        Fv         = scatteredInterpolant(P(:,1)*f , ... % N-coordinate
                     P(:,2), ...                         % Z-coordinate
                     S, 'natural', 'none');              % Salinity
        
        Varv = Fv(nq*f , zq);                            % Interpolate as mesh-coordinates
        
        Fv_ex      = scatteredInterpolant(P(:,1)*f , ... % N-coordinate
                     P(:,2), ...                         % Z-coordinate
                     S, 'natural', 'nearest');              % Salinity    
        Varv_ex = Fv_ex(nq*f , zq);   
        
        msh{ct}.CTD_int.Sal(:,:,hour,1) = Varv;   
        msh{ct}.CTD_int.Sal_ex(:,:,hour,1) = Varv_ex;
        msh{ct}.CTD_int.Z(:,:,hour,1)   = zq;            
        msh{ct}.CTD_int.X(:,:,hour,1)   = nq;            
        
    end
end

% save('MSH_20150914', 'MSH_20150914')
%% Figure of salinity in cross-sections
msh = MSH_20150914;
% Create colormap
cmap = brewermap(15, 'YlOrBr'); % Brewermap.m (available through Github, use the same code); map = brewermap(N,scheme)
colormap(cmap)     

plotlist = 1:1:13;

for ct=1:3
    
    figure(ct)
    
    for hour=1:13

        subplot(7,2,plotlist(hour)); 

        patch(msh{ct}.p.N(:,:,1),...    % Plot salinity for N, Z
              msh{ct}.p.Z(:,:,1),...
              msh{ct}.CTD_int.Sal_ex((msh{ct}.p.progfgood(:,hour)))', 'LineStyle', 'none', 'FaceColor', 'flat'); colormap(cmap); colorbar           

        caxis([0,25]); hold on
        plot((msh{ct}.p.nbed),...       % Plot channel bed
              msh{ct}.p.zbed,'-k', 'LineWidth', 1.15);
        hold off
        
    end
end


% %% Example figure: x-velocity
% 
% %tid_adcp = [1:13;    % Transect used to create a figure (for the second transect, there are some partial transect which are left out)
%  %   1,2,3,6:15;
%   %  1:13];
%   
% tid_adcp = [1:13;    % Transect used to create a figure (for the second transect, there are some partial transect which are left out)
%      1:13;
%     1:13];
%   
% dir           = [-1, 1, 1]; % Direction of flow
% 
% position      = [1,0, 22, 30];        % Size and position of figure
% no_figs1      = 6;                    % Total number of figures (no_figs1 times no_figs2)
% no_figs2      = 6;
% ax_v          = 1.5;                  % Size of subplot: height
% ax_h          = 4;                    % Size of subplot: width
% s             = 0.2;
% 
% 
% fig           = figure;
% u            = fig.Units;
% fig.Units    = 'centimeters';
% set(fig, 'Position', position);
% % set(gcf, 'Position', get(0, 'Screensize'));
% cmap = brewermap(20, 'RdBu');
% %%%%%%%%%%%%%%%%%%% Plots with Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for ct = 1:3
%     for t = 1:13
%         fig          = subplot(no_figs1,no_figs2,36);
%         fig.Units    = 'centimeters';
%         set(fig, 'Position', [2+(ct-1)*(s+ax_h), position(4)-(1+t)*1*(s+ax_v)-3, ax_h, ax_v]);
%         
%         %%%%%%%%%%%%% Use function 'patch' to plot the grid with velocities
%         p        = patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),...
%             dir(ct)* msh{ct}.pars(msh{ct}.p.progfgood_vec(:,tid_adcp(ct,t),1))');
%         %%%%%%%%%%%%%
%         hold on % This enable plotting something else 'on top of' the first plot.
%         plot(msh{ct}.p.nbed, msh{ct}.p.zbed, '-k')% add bed geometry
%         caxis([-1,1]);              % Set the limits of velocity-values
%         set(p, 'edgecolor','none')  % Removes black lines of grid
%         ylim([-15 2])
%         xlim([-1 350])
%         colormap(cmap)
%         if ct==1 % Here, t=# is added for each row
%             handle=title(['\fontsize{8} t = ', num2str(t)],'FontWeight','Normal'); handle.Units='centimeters'; set(handle,'Position',[-1.5, 0.5, 1]);
%             ylabel('Depth [m]')
%         else
%             set(gca,'YTickLabel',[]);
%         end
%         
%         if t~=13
%             set(gca,'XTickLabel',[]);
%         else
%             xlabel('x-coordinate [m]')
%         end
%         %axis tight
%         
%         if t==1 && ct==2
%             handle=title(['\fontsize{10} Flow velocity in x-direction'],'FontWeight','Normal'); handle.Units='centimeters'; set(handle,'Position',[1.3,1.6, 1]);
%         end
%     end
% end
% % Include legend:
% lgd       = colorbar;  
% lgd.Units = 'centimeters'; 
% set(lgd, 'Position', [2.2+(3)*(s+ax_h), position(4)-(3)*(s+ax_v)-3, 0.4, 2*ax_v]);
% lgd.Label.String = 'u [m s^{-1}]' ;
% % 
% % %export_fig([ROOTFOLDER, '\Figures\Velocity'], '-png', '-transparent', '-r464')
% 
