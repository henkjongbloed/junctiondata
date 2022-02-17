%% set preferences
nr_transects    = 3;     % [#] Number of transect in the dataset. 
mesh_size_x     = 10;    % [m] Size of the grid: width of columns; can be changed
mesh_size_z     = 1.5;   % [m] Size of the grid: height of rows; can be changed

% load data - water levels measured at station Spijkenisse (downloaded from
% https://waterinfo.rws.nl/#!/nav/index/)
% Time in UTC+1 (MET)

load 20140812_h_Spijkenisse.mat % h.date h.level 

%% read in and filter raw adcp data (.000, .dat, .txt)
% 
adcp            = cell(3,1);                                    % Initialize empty cell array. In each cell the data from one transect will be stored. 

% ADCP data is read through function 'readDeployment'  
DEPNAME = 'OMHA'; % readDeployments read files starting with DEPNAME

%load('adcp600.mat')
for ct = 1:3
PATH    = [ROOTFOLDER,'\Data\master',...
    '\raai',num2str(ct), '\transect'];                          % reads files in PATH. PATH wordt gemaakt door een pad naar de folder waar de data staat 
                                                                %(zou bij jou moeten eindigen in Data\Aquavision_ADCP_master_600. 
                                                                %Vervolgens wordt '\raai + cijfer + \transect toegevoegd                                                 
adcp{ct}        = readDeployment(DEPNAME,PATH); % Nu wordt alle data ingelezen, komt in workspace terecht.
end

%% Interpolate data every hour to mesh 'msh'

msh = cell(3,1); % Initialize empty cell array to store the results
        
for ct = 1:3
    time    = datenum(adcp{ct}.timeV);
    eta     = interp1(h.date,h.level,time); % Calculate the waterlevel at each time step (using interpolation)
    msh{ct} = procTrans(...
        adcp{ct},adcp{ct}.FileNumber,'DepthTransducer',0.3,'DeltaN', mesh_size_x, 'DeltaZ', mesh_size_z, 'Eta',eta,...
        'MinimumSigma',0.06,'CumulateCrossings',false,'ConventionalProcessing',false,'TopMeshLowestEta', true,...
        'ConstantZetaMesh', true, 'RemoveOutliers', 0, 'Proximity', 0, 'StdFiltering', 6,'ShipReference','bt',...
        'RotatePars', true, 'EnableDebugging',false); 
end

% save('WURData_processed.mat', 'msh', 'adcp') % Can save the results for later use. 
