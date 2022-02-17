%% Read, filter, plot ADCP data
% Created by Henk Jongbloed

clearvars % Clear all variables from the workspace
close all % Close all figures
clc
RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'\DataAnalysis'])
DS = 'NMOMNW15'; %DataSet: 'OMHA14' or 'NMOMNW15'
addPaths(RF, DS, 'old')
%% Import preprocessed semi-raw data (ADCP / CTD / H structs)
load(['ADCP',DS,'.mat']);
h = load(['H',DS,'.mat']);
h.time            = datetime(h.date, 'ConvertFrom','datenum','Format','HH:mm:ss');
if strcmp(DS,'NMOMNW15')
    ctd = load(['CTD',DS,'.mat']).CTD_20150914 ;
end
eta = h.level; date0 = h.date(1); lat = 51.92; interval = 24*(h.date(2)-h.date(1));
[tidestruct,etaout] = t_tide(eta,...
    'interval',interval, ...        % hourly data
    'start',date0,...               % start time is datestr(tuk_time(1))
    'latitude',lat,...              % Latitude of obs
    'error','linear',...            % coloured boostrap CI
    'synthesis',1);                 % Use SNR=1 for synthesis.


ctd = hctd(h,ctd,DS);

%% Interpolate data every hour to mesh 'msh'
transects    = 1;     % [#] Number of transect in the dataset.
mesh_size_x     = 15;    % [m] Size of the grid: width of columns; can be changed
mesh_size_z     = .75;   % [m] Size of the grid: height of rows; can be changed
msh = cell(numel(transects),1); % Initialize empty cell array to store the results
%load('C:\Users\jongb013\Documents\PHD\2-Programming\DataAnalysis\processedData\H_2015091415.mat');
%h=H_2015091415;

%% Process
% figure(1)
% for ct = 1:nr_transects
%     [x,y] = utmADCP(adcp{ct});
%     hold on
%     plot(x,y)
% end

for ct = transects
    time    = datenum(adcp{ct}.timeV);
    eta     = interp1(h.date,h.level, time); % Calculate the waterlevel at each time step (using interpolation)
    msh{ct} = procTrans(...
        adcp{ct},adcp{ct}.FileNumber,'ConventionalProcessing',false,...
        'ConstantZetaMesh',true, 'TopMeshLowestEta', true,...
        'DeltaN', mesh_size_x, 'DeltaZ', mesh_size_z, 'Eta', eta,...
        'CumulateCrossings',false,'ShipReference','bt','EnableDebugging',false);
    msh{ct}.adcp_good = 1:13;
end

%% Interpolate CTD Data
%[h,ctd] = hctd(h,ctd,DS);

for ct = transects
    % start hour (time) loop
    for hour = 1:13
        %         X=[];Y=[];T=[];Z=[];S=[];    % Initiate variables
        
        % for every hour and transect find indices
        ind        = find(and(ctd.Hour(:) == hour, ctd.Trans(:) == ct)); %Aangepast
        X          =  ctd.X(ind);
        Y          =  ctd.Y(ind);
        Z          =  ctd.Z_nap(ind);
        %Z          =  ctd.Z(ind);
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


slider_plot(msh)
%% Velocity: QUIVER
%close all
% msh = MSH_20140812;
% dir           = [-1, 1, -1]; % Direction of flow
% velmap = brewermap(20, 'RdBu');
% velFig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
% i=1;
% hourV = [1:13];
% ctV = [1 2 3];
% %hourV = [12];
% %ctV = [2];
% clear mav miv mav2 miv2
%
% mav = 0; miv = 0;
% for hour = hourV
%     for ct=ctV
%         vel = dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood(:,msh{ct}.adcp_good(hour)))';
%
%         mav2 = max(max(vel)); miv2 = min(min(vel));
%         mav = max(mav,mav2);
%         miv = min(miv,miv2);
%         disp(['At t = ',num2str(hour), ' we have max = ', num2str(mav2),' and min = ', num2str(miv2)])
%         i=i+1;
%     end
% end
%
%
%
%
% i=1;
%
% for hour = hourV
%     figure(hour);
%     for ct = ctV
%         titles = {'NM', 'OM', 'RWW'};
%         subplot(length(ctV), 1, ct)
%         title([titles{ct},': hour = ', num2str(hour)])
%
%         %caxis([-1,1]);
%         hold on
%         caxis([-1.5,1.5]);
%         vel = dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood(:,msh{ct}.adcp_good(hour)))';
%
%         patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),vel, 'LineStyle', '-', 'LineWidth',.2, 'FaceColor', 'flat');
%
%         quiver(msh{ct}.N, msh{ct}.Z(:,:,1),...
%             msh{ct}.sec.pars(:,:,hour,2),...
%             msh{ct}.sec.pars(:,:,hour,3), 'color', 'k', 'LineWidth',1.5,'ShowArrowHead','on','AutoScale','on');
%         plot((msh{ct}.p.nbed),...       % Plot channel bed
%             msh{ct}.p.zbed,'-k', 'LineWidth', 1.15); colormap(velmap); cb = colorbar;
%         hold off
%         %         caxis([0,25]);
%         % axis equal
%         axis tight
%         %         print(sp)
%         i=i+1;
%     end
% end
%
%
%
% % for hour = hourV
% %     for ct = ctV
% %         %hour = 8;
% %         %ct = 2;
% %         titles = {'NM', 'OM', 'RWW'};
% %         %subplot(13,2, 2*(hour-1)+ct)
% %         if length(ctV) == 1
% %             if any(length(hourV) == [1 2])
% %                 subplot(length(hourV),length(ctV),i)
% %             elseif any(length(hourV) == [3, 4])
% %                 subplot(2,2,i)
% %             elseif any(length(hourV) == [5, 6])
% %                 subplot(2,3,i)
% %             elseif any(length(hourV) == [7, 8, 9])
% %                 subplot(3,3,i)
% %             elseif length(hourV) == 13
% %                 subplot(3,5,i)
% %             else
% %                 subplot(3,4,i)
% %             end
% %         else
% %             subplot(length(hourV),length(ctV),i)
% %         end
% %         title([titles{ct},': hour = ', num2str(hour)])
% %
% %         %caxis([-1,1]);
% %         hold on
% %         caxis([-1.5,1.5]);
% %         vel = dir(ct)* msh{ct}.sec.pars(msh{ct}.p.progfgood(:,msh{ct}.adcp_good(hour)))';
% %
% %         patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1),vel, 'LineStyle', '-', 'LineWidth',.2, 'FaceColor', 'flat');
% %
% %         quiver(msh{ct}.N, msh{ct}.Z(:,:,1),...
% %             msh{ct}.sec.pars(:,:,hour,2),...
% %             msh{ct}.sec.pars(:,:,hour,3), 'color', 'k', 'LineWidth',1.5,'ShowArrowHead','on','AutoScale','on');
% %         plot((msh{ct}.p.nbed),...       % Plot channel bed
% %             msh{ct}.p.zbed,'-k', 'LineWidth', 1.15); colormap(velmap); cb = colorbar;
% %         hold off
% %         %         caxis([0,25]);
% %         % axis equal
% %         axis tight
% %         %         print(sp)
% %         i=i+1;
% %     end
% % end
%
%
%
% sgtitle('Velocity, rotated to match with the direction of the cross-section')
%
%
%
% %% Salinity
%
% salmap = brewermap(15, 'YlOrBr'); % Brewermap.m (available through Github, use the same code); map = brewermap(N,scheme)
% %colormap(salmap)
%
% salFig = figure('units','normalized','outerposition',[0 0 1 1]);
% i=1;
%
% mas = 0; mis = 0;
% for hour=hourV
%     for ct=ctV
%         sal = msh{ct}.CTD_int.Sal_ex((msh{ct}.p.progfgood(:,hour)))';
%         mas = max(mas,max(max(sal)));
%         mis = min(mis,min(min(sal)));
%         i = i+1;
%     end
% end
%
% i=1;
%
% for hour=hourV
%     figure(hour);
%     for ct=ctV
%         titles = {'NM', 'OM', 'RWW'};
%         subplot(length(ctV), 1, ct)
%         title([titles{ct},': hour = ', num2str(hour)])
%         sal = msh{ct}.CTD_int.Sal_ex((msh{ct}.p.progfgood(:,hour)))';
%         patch(msh{ct}.p.N(:,:,1),...    % Plot salinity for N, Z
%             msh{ct}.p.Z(:,:,1),sal, 'LineStyle', 'none','FaceColor', 'flat'); colormap(salmap); cb = colorbar;
%
%         caxis([mis, mas]);
%         hold on
%         plot((msh{ct}.p.nbed),...       % Plot channel bed
%             msh{ct}.p.zbed,'-k', 'LineWidth', 1.15);
%         hold off
%         i = i+1;
%     end
% end
%
%
%
%
% % for hour=hourV
% %     figure;
% %     for ct=ctV
% %         if length(ctV) == 1
% %             if any(length(hourV) == [1 2])
% %                 subplot(length(hourV),length(ctV),i)
% %             elseif any(length(hourV) == [3, 4])
% %                 subplot(2,2,i)
% %             elseif any(length(hourV) == [5, 6])
% %                 subplot(2,3,i)
% %             elseif any(length(hourV) == [7, 8, 9])
% %                 subplot(3,3,i)
% %             elseif length(hourV) == 13
% %                 subplot(3,5,i)
% %             else
% %                 subplot(3,4,i)
% %             end
% %         else
% %             subplot(length(hourV),length(ctV),i)
% %         end
% %         title([titles{ct},': hour = ', num2str(hour)])
% %         sal = msh{ct}.CTD_int.Sal_ex((msh{ct}.p.progfgood(:,hour)))';
% %         patch(msh{ct}.p.N(:,:,1),...    % Plot salinity for N, Z
% %             msh{ct}.p.Z(:,:,1),sal, 'LineStyle', 'none','FaceColor', 'flat'); colormap(salmap); cb = colorbar;
% %
% %         caxis([mis, mas]);
% %         hold on
% %         plot((msh{ct}.p.nbed),...       % Plot channel bed
% %             msh{ct}.p.zbed,'-k', 'LineWidth', 1.15);
% %         hold off
% %         i = i+1;
% %     end
% % end
%
%
%
%
%
%
% %colormap(salmap);
% %cb = colorbar;
% % cb = colorbar; cb.Units = 'centimeters'; %('Position',[
% %% Analysis of transport (does not work anymore)
% % msh = compQ_v2(msh);
% % for ct=1:nr_transects
% %     msh{ct}.Qloc = zeros(size(squeeze(msh{ct}.pars(:,:,:,1))));
% %     msh{ct}.Q = nan(size(squeeze(msh{ct}.pars(1,1,:,1))));
% %     msh{ct}.adcp_good = 1:13;
% %     for hour=1:13
% %         %x,y,z coordinates
% %         vdotn = squeeze(msh{ct}.pars(:,:,msh{ct}.adcp_good(hour),1))  * msh{ct}.Nvec(1) +...
% %             squeeze(msh{ct}.pars(:,:,msh{ct}.adcp_good(hour),2))  * msh{ct}.Nvec(2) +...
% %             squeeze(msh{ct}.pars(:,:,msh{ct}.adcp_good(hour),3))  * 0;
% %         Qloc = msh{ct}.A(:,:,1).*vdotn;
% %         msh{ct}.Qloc(:,:,msh{ct}.adcp_good(hour)) = Qloc;
% %         msh{ct}.Q(msh{ct}.adcp_good(hour)) = nansum(Qloc,[1,2]);
% %         %          msh{ct}.Qloc(:,:,msh{ct}.adcp_good(hour)) = Qloc;
% %         Qlocrot = msh{ct}.sec.pars(:,:,msh{ct}.adcp_good(hour),1).*msh{ct}.A(:,:,1);
% %         msh{ct}.sec.Qloc(:,:,msh{ct}.adcp_good(hour)) = Qlocrot;
% %         msh{ct}.sec.Q(msh{ct}.adcp_good(hour)) = nansum(Qlocrot,[1,2]);
% %         Qlocrot_S = Qlocrot.*msh{ct}.CTD_int.Salinity_ex(:,:,hour);
% %         msh{ct}.sec.QlocS(:,:,msh{ct}.adcp_good(hour)) = Qlocrot_S;
% %         msh{ct}.sec.QS(msh{ct}.adcp_good(hour)) = nansum(Qlocrot_S,[1,2]);
% %     end
% % end%% bin
% % year = '2015';
% % [h, msh] = importHMSH(year); %imports a pre-processed structure MSH with (ReadDeployment and ProcTrans already applied)
% %
% % %% Tidal Analysis
%
%
% %% Possible operations on data
% %msh = comp_plotQ(msh);
% %% read in and filter raw adcp data (.000, .dat, .txt)
%
% % adcp            = cell(3,1);                                    % Initialize empty cell array. In each cell the data from one transect will be stored.
% %
% % % ADCP data is read through function 'readDeployment'
% % DEPNAME = 'OSR'; % readDeployments read files starting with DEPNAME
% %
% % %load('adcp600.mat')
% % for ct = 1:3
% % PATH    = [RF, 'Data\', DS,'\data_meting_14sep2015',...
% %     '\raai',num2str(ct)];                          % reads files in PATH. PATH wordt gemaakt door een pad naar de folder waar de data staat
% %                                                                 %(zou bij jou moeten eindigen in Data\Aquavision_ADCP_master_600.
% %                                                                 %Vervolgens wordt '\raai + cijfer + \transect toegevoegd
% % adcp{ct} = readDeployment(DEPNAME,PATH); % Nu wordt alle data ingelezen, komt in workspace terecht.
% % end