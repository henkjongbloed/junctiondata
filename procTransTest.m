%procTransTest 15-9-20

clearvars % Clear all variables from the workspace
close all % Close all figures
clc
%% set preferences
%cd('C:\Users\jongb013\OneDrive - WageningenUR\PhDHenkJongbloed\2. Programming\Code\DataAnalysis')
addPaths

%load H_2015091415.mat % Imports waterlevels Geulhaven, 14-Sep-15 00:00 until 15-Sep-15 12:24, 1 min interval
load H_201408.mat
%load CTD_20150914.mat  % Pre-processed semi-raw data from CTD sensor
load ADCP_20140812.mat
% load MSH_20150914.mat %1: NM 2: OM 3: RWW (salt and velocity data)

h = H_201408;
adcp = ADCP_20140812;
%% Interpolate data every hour to mesh 'msh'

msh = cell(1,1); % Initialize empty cell array to store the results

%nr_transects    = 3;     % [#] Number of transect in the dataset. -> not used
mesh_size_x     = 3;    % [m] Size of the grid: width of columns; can be changed
mesh_size_z     = 1;   % [m] Size of the grid: height of rows; can be changed

for ct = 1:3
    time    = datenum(adcp{ct}.timeV);
    eta     = interp1(h.date,h.level,time); % Calculate the waterlevel at each time step (using interpolation)
    msh{ct} = procTrans(...
        adcp{ct},adcp{ct}.FileNumber,'DepthTransducer',0.3,'DeltaN', mesh_size_x, 'DeltaZ', mesh_size_z, 'Eta',eta,...
        'MinimumSigma',0.06,'CumulateCrossings',false,'ConventionalProcessing', false, 'TopMeshLowestEta', true,...
        'ConstantZetaMesh', true, 'RemoveOutliers', 0, 'Proximity', 0, 'StdFiltering', 6,'ShipReference','bt',...
        'RotatePars', true, 'EnableDebugging',false);
    disp(ct)
end
save('MSH_2014Test', 'msh')


%%
load MSH_20150914.mat 
load H_2015091415.mat
load CTD_20150914.mat
                                                    % only time values
H_2015091415.time            = datetime(H_2015091415.date, 'ConvertFrom','datenum','Format','HH:mm:ss');
                                                    % Change time CTD data
                                                    % from UTM --> UTM+1
hour              = datenum('01-Jan-2015 01:00:00')-datenum('01-Jan-2015 00:00:00'); %hour*24 = 1
CTD_20150914.time = datetime(CTD_20150914.Time2+hour, 'ConvertFrom','datenum','Format','HH:mm:ss');
H_2015091415.datestr         = H_2015091415.datestr(301:1440);              % Create subset h-data, so that time only occurs once and corresponds to time of CTD-data. 
H_2015091415.time            = H_2015091415.time(301:1440);
H_2015091415.date            = H_2015091415.date(301:1440);
H_2015091415.level           = H_2015091415.level(301:1440);
% Divide data per transect and hour and vertical... 
%  Interpolation waterlevels Geulhaven
CTD_20150914.time_unique = unique(CTD_20150914.time);
H_2015091415.level_int     = interp1(H_2015091415.time, H_2015091415.level, CTD_20150914.time_unique, 'pchip'); %Piecewise Cubic Hermite Interpolating Polynomial.
H_2015091415.time_int      = CTD_20150914.time_unique;
CTD_20150914.h          = interp1(H_2015091415.time, H_2015091415.level, CTD_20150914.time, 'pchip'); %Piecewise Cubic Hermite Interpolating Polynomial.
CTD_20150914.Z_nap      = -CTD_20150914.Z + CTD_20150914.h;
%% Interpolate salt profiles per hour
%figure;
% msh = MSH_20150914;
% Start transect loop
for ct =1:3
    % start hour (time) loop
    for hour = 1:13
%         X=[];Y=[];T=[];Z=[];S=[];    % Initiate variables
        
        % for every hour and transect find indices
        ind        = find(and(CTD_20150914.Hour(:) == hour, CTD_20150914.Trans(:) == ct)); %Aangepast
        X          =  CTD_20150914.X(ind);
        Y          =  CTD_20150914.Y(ind);
        Z          =  CTD_20150914.Z_nap(ind);
        S          =  CTD_20150914.S(ind);
       
        %%%%%%%%% Step 1.  Convert transect to to n,z-coordinates %%%%%%%%%
        %%%%%%%%% -> what about s?
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

%INCLUDE: Understand Model and Known etc...

%         'DepthTransducer'
%           scalar value - Sets the depth of the transducers in m. Default
%           is 0.3m.
%         'DeltaN'
%           scalar value - Sets the mesh size in n-direction. Default is 5m
%         'DeltaZ'
%           scalar value - Sets the mesh size in z-direction. Default is 1m
%         'Eta'
%           1xE row vector - Gives the value of eta (Water level) for each
%           ensmeble. Default is 0 (no water level changes)
%         'MinimumSigma'
%           scalar value between 0 and 1 - Sets the minium sigma to start
%           the mesh. Default is 0.06.
%         'CumulateCrossings'
%           true | {false} - Whether the time-steps given in TID for each
%           row should be cumulated
%         'ConventionalProcessing'
%           true | {false} - Whether to use the conventional processing.
%         'TopMeshLowestEta'
%           true | {false} - Whether to set the top of the mesh to the
%           lowest value of eta or to the top of the data thoughout the
%           dataset
%         'ConstantZetaMesh'
%           true | {false} - Whether to use a mesh constant in Z
%           coordinates instead of constant in Sigma coordinates
%         'RemoveOutliers'
%           positive scalar - Indicates how many times the residual in
%           beam-velocity should exceed the median of all residuals to 
%           discard a beam-velocity from the fit. Default is 0, indicating 
%           no outlier removal
%         'Proximity'
%           positive scalar - How close should the beam-velocity be in 
%           s-direction (i.e. along Nvec) to be included in the fit. 
%           Default is 0, indicating to include all data
%         'StdFiltering'
%           positive scalar - Indicates how many times the standard
%           deviation of any parameter in a cell should excess the median
%           of the standard deviations over all cells to be treated as bad.
%           Default is 6. A value of 0 deactivates this filtering.
%         'ShipReference'
%           {'bt'} | 'gps' | 'btgps' | 'gpsbt' - Indicates which method
%           should be used to compute the ship velocity:
%           'bt' use bottom-tracking
%           'gps' use GPS
%           'btgps' use bottom-tracking and GPS if bottom-tracking is
%               unavailable
%           'gpsbt' use GPS and bottom-tracking when GPS is unavailable
%         'Pusr'
%           two element row vector - Sets a user defined center of the
%           cross-section. Can be usefull when the measurements are not
%           uniformly collected along the cross-section
%         'RotatePars'
%           {true} | false - Whether to perform rotations
%         'Model'
%           scalar function handle - provides a function to specify a
%           velocity model. The function should accept as input nx5 matrix
%           where each column represents the n, z, sigma, s and t 
%           coordinates of a beam velocity. The function must return three
%           matrices representing the models for u, v and w respectively.
%           The matrices have n rows and as many columns as parameters in
%           the model. Default function fits a constant mean velocity for
%           u, v and w at each time step.
%         'Known'
%           scalar function handle - provides a function to specify the
%           know terms in the velocity model. The function should accept as
%           input nx5 matrix where each column represents the n, z, sigma, 
%           s and t coordinates of a beam velocity. The function must 
%           return three vectors representing the known terms for u, v and 
%           w respectively. 
%         'GetVelocity'
%           scalar function handle - defines how to retrieve velocity from
%           the estimated parameters. It accepts and nxP matrix as input.
%           It returns three matrices with the estimated u,v and w
%           velocity.

%         Instead of defining the velocity models through function handles
%         it is also possible to define them with the following parameters
%         (included for backward compatibility):
%         'ModelU_t'
%         'ModelV_t'
%         'ModelW_t'
%           PxE matrix - Each row corresponds to a model parameter for the
%           velocity component. The values are the time-varying 
%           coefficients of the velocity model. Default is ones(1,E) which
%           corresponds to computing the average velocity.
%         'KnownU_t'
%         'KnownV_t'
%         'KnownW_t'
%           Px1 row-vector - Allows to set any known term of the velocity 
%           model. Default is zeros(1,E) which corresponds to no known term 
%         Note that if this last way of defining a velocity model is used
%         the functions given in the 'Model', 'Known' and
%         'GetVelocity' inputs will be ignored.


% msh = procTrans(adcp,tid,varargin);