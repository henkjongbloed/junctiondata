%% Read, filter, plot ADCP data 
% for Judith Poelman HbR, Wur
% 2017-04-11 MC Verbeek

%% example to read in ADCP data 14 sep

% SET PATHS: [SET YOUR OWN]
% Path_ADCPtools  = 'D:\_data\adcptools_151112';
% Path_ADCPdata   = 'D:\mverbeek2\surfdrive_files\_administratie\partners\Poelman_Judith_HbR\data_meting_14sep2015\';
% 
% addpath         Path_ADCPtools
% addpath         Path_ADCPdata
RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
DS = 'NMOMNW15'; %DataSet: 'OMHA14' or 'NMOMNW15'
addPaths(RF,DS)
% set preferences
nr_transects    = 3;     % [#]
mesh_size_x     = 20;    % [m]
mesh_size_z     = 1.5;   % [m]

% load data - water level 'eta' geulhaven in MET (=UTC-1)
load            20150914_h_geulhaven;   

% read in and filter raw adcp data (.000, .dat, .txt)
% adcp            = cell(3,1);
% for ct = 1:3
% adcp{ct}        = readDeployment('OSR',[Path_ADCPdata,'\raai',num2str(ct)]);
% end
load ADCPNMOMNW15.mat
%% interpolate data every hour to mesh 'msh'
hour = datenum('01-Jan-2015 01:00:00')-datenum('01-Jan-2015 00:00:00');
for ct = 1:nr_transects
    time    = datenum(adcp{ct}.timeV)-repmat(hour,[length(adcp{ct}.timeV),1]);
    eta     = interp1(h.date,h.level,time);
    msh{ct,1}=procTrans(...
        adcp{ct},adcp{ct}.FileNumber,'ConventionalProcessing',false,...
        'ConstantZetaMesh',true, 'TopMeshLowestEta', true,...
        'DeltaN', mesh_size_x, 'DeltaZ', mesh_size_z, 'Eta',eta,...
        'CumulateCrossings',false,'ShipReference','bt','EnableDebugging',false);
end

%% make plots adcp data per channel

printje     = 0;        % [0 or 1], if 1 plots are saved

% plot velocity data channel 1

ct          = 1;        % channel number
naam        = '20150930_transect1_adcp.pdf';
Ts_1        = 1:14;     % time steps (hours)
figure(1)
for j = Ts_1
subaxis(7,2,find(Ts_1==j), 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
h = patch(msh{ct}.p.N(:,:,1), msh{ct}.p.Z(:,:,1), msh{ct}.cs.pars(msh{ct}.p.progfgood(:,j))'); colorbar
hold on

fn = 50; fz = 50;
p = quiver(msh{ct}.N(msh{ct}.p.fgood),...  % channel-normal positions
           msh{ct}.Z(msh{ct}.p.fgood),...  % channel-vertical positions
           msh{ct}.cs.pars(msh{ct}.p.progfgood_vec(:,j,2)).*fn,...
           msh{ct}.cs.pars(msh{ct}.p.progfgood_vec(:,j,3)).*fz,...
           'DisplayName', 'secondary flow',...
           'AutoScale','off', 'color',[0 0 0], ...
           'MaxHeadSize', 0);

plot(msh{ct}.p.nbed, msh{ct}.p.zbed, '-k')% add bed geometry
caxis([-1.5,1.5])
text(msh{ct}.p.nbed(1), msh{ct}.p.zbed(1)+1.5, ['t= ', num2str(round((msh{ct}.time(j)-msh{ct}.time(1))*24+1))])
    axis tight
    axis off
set(h, 'edgecolor','none')
end
hold off
set(gcf,'paperpositionmode','auto');
if printje
print(gcf, '-dpdf', '-r600', naam);
end

% plot velocity data channel 2

figure(2)
ct          = 2;
Ts_1        = 1:14;
naam        = '20150930_transect2_adcp.pdf';

for j = Ts_1
subaxis(7,2,find(Ts_1==j), 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
h= patch(msh{ct,1}.p.N(:,:,1), msh{ct,1}.p.Z(:,:,1), msh{ct,1}.cs.pars(msh{ct,1}.p.progfgood(:,j))'); colorbar
hold on
plot(msh{ct,1}.p.nbed, msh{ct,1}.p.zbed, '-k')
text(msh{ct,1}.p.nbed(1), msh{ct,1}.p.zbed(1)+1.5, ['t= ', num2str(round((msh{ct,1}.time(j)-msh{ct,1}.time(1))*24+1))])
    axis tight
    axis off
caxis([-1.5,1.5])
set(h, 'edgecolor','none')
end
hold off

set(gcf,'paperpositionmode','auto');
if printje
print(gcf, '-dpdf', '-r600', naam);
end

% plot velocity data channel 3

figure(3)
ct      = 3;
Ts_1    = [1:9,11:13,15];
naam    = '20150930_transect3_adcp.pdf';

for j = Ts_1
subaxis(7,2,find(Ts_1==j), 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
h= patch(msh{ct,1}.p.N(:,:,1), msh{ct,1}.p.Z(:,:,1), msh{ct,1}.cs.pars(msh{ct,1}.p.progfgood(:,j))'); colorbar
hold on
plot(msh{ct,1}.p.nbed, msh{ct,1}.p.zbed, '-k')
text(msh{ct,1}.p.nbed(1), msh{ct,1}.p.zbed(1)+1.5, ['t= ', num2str(round((msh{ct,1}.time(j)-msh{ct,1}.time(1))*24+1))])
    axis tight
    axis off
caxis([-1.5,1.5])
set(h, 'edgecolor','none')
end
hold off

set(gcf,'paperpositionmode','auto');
if printje
print(gcf, '-dpdf', '-r600', naam);
end

