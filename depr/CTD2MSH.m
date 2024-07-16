clearvars
close all
clc
addPaths
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
msh = MSH_20150914;
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
        
        Salinity          =  CTD_20150914.S(ind);
        Temperature          =  CTD_20150914.T(ind);
        Pressure   = CTD_20150914.Pressure(ind);
        SC         = CTD_20150914.SC(ind);
        SoundV     = CTD_20150914.SoundV(ind);
        Density    = CTD_20150914.Density(ind);
        Z_nap      = CTD_20150914.Z_nap(ind)  ;
        %%%%%%%%% Step 1.  Convert transect to to n,z-coordinates %%%%%%%%%
        % n - along the transect, z - upwards
        [x,y]      = wgs2utm(Y,X);
        
        P_tmp      = [x,y,Z];                            % create CTD position data vector
         
        % Transformation of the coordinates to the mesh; following the method that was applied to the ADCP data (in procTrans2.m by B. Vermeulen) 
        P1          = msh{ct}.Tvec(1)*(P_tmp(:,1)-msh{ct}.Pm(1))+msh{ct}.Tvec(2)*(P_tmp(:,2)-msh{ct}.Pm(2));
        P          = [P1, P_tmp(:,3)];

        %%%%%%%%% Step 2. Interpolate data using scattered interpolant %%%%
        % Define mesh location at which interpolant will be evaluated

        nq         = msh{ct}.N(:,:,1)-max(max(msh{ct}.N(:,:,1)))./2; % N values in msh, transformed as Pm should be centre of the mesh. 
        zq         = msh{ct}.Z(:,:,1);
        
        % Define and plot interpolant ------------------------------------------
        f          = 1/50;%1/30;% Why 30?
        % Salt
        Fv         = scatteredInterpolant(P(:,1)*f , ... % N-coordinate
                     P(:,2), ...                         % Z-coordinate
                     Salinity, 'natural', 'none');              % Salinity
        
        Varv = Fv(nq*f , zq);                            % Interpolate as mesh-coordinates
        
        Fv_ex      = scatteredInterpolant(P(:,1)*f , ... % N-coordinate
                     P(:,2), ...                         % Z-coordinate
                     Salinity, 'natural', 'nearest');              % Salinity    
        Varv_ex = Fv_ex(nq*f , zq);   
        
        msh{ct}.CTD_int.Salinity(:,:,hour,1) = Varv;   
        msh{ct}.CTD_int.Salinity_ex(:,:,hour,1) = Varv_ex;
        msh{ct}.CTD_int.Z(:,:,hour,1)   = zq;            
        msh{ct}.CTD_int.X(:,:,hour,1)   = nq;
        
        % Temperature
        Fv.Values = Temperature   ;
        Varv = Fv(nq*f , zq);                            % Interpolate as mesh-coordinates
        Fv_ex.Values = Temperature; 
        Varv_ex = Fv_ex(nq*f , zq); 
        
        msh{ct}.CTD_int.Temperature(:,:,hour,1) = Varv;   
        msh{ct}.CTD_int.Temperature_ex(:,:,hour,1) = Varv_ex;
        
        % Density
        Fv.Values = Density   ;
        Varv = Fv(nq*f , zq);                            % Interpolate as mesh-coordinates
        Fv_ex.Values = Density; 
        Varv_ex = Fv_ex(nq*f , zq);
        msh{ct}.CTD_int.Density(:,:,hour,1) = Varv;   
        msh{ct}.CTD_int.Density_ex(:,:,hour,1) = Varv_ex;
        
        % SoundV
        Fv.Values = SoundV   ;
        Varv = Fv(nq*f , zq);                            % Interpolate as mesh-coordinates
        Fv_ex.Values = SoundV; 
        Varv_ex = Fv_ex(nq*f , zq);
        msh{ct}.CTD_int.SoundV(:,:,hour,1) = Varv;   
        msh{ct}.CTD_int.SoundV_ex(:,:,hour,1) = Varv_ex;        
          
        % Pressure
        Fv.Values = Pressure   ;
        Varv = Fv(nq*f , zq);                            % Interpolate as mesh-coordinates
        Fv_ex.Values = Pressure; 
        Varv_ex = Fv_ex(nq*f , zq);
        msh{ct}.CTD_int.Pressure(:,:,hour,1) = Varv;   
        msh{ct}.CTD_int.Pressure_ex(:,:,hour,1) = Varv_ex;
        
        % SC -> what is this??
        Fv.Values = SC   ;
        Varv = Fv(nq*f , zq);                            % Interpolate as mesh-coordinates
        Fv_ex.Values = SC; 
        Varv_ex = Fv_ex(nq*f , zq);
        msh{ct}.CTD_int.SC(:,:,hour,1) = Varv;   
        msh{ct}.CTD_int.SC_ex(:,:,hour,1) = Varv_ex;
        
        % Z_nap
        Fv.Values = Z_nap   ;
        Varv = Fv(nq*f , zq);                            % Interpolate as mesh-coordinates
        Fv_ex.Values = Z_nap; 
        Varv_ex = Fv_ex(nq*f , zq);
        msh{ct}.CTD_int.Z_nap(:,:,hour,1) = Varv;   
        msh{ct}.CTD_int.Z_nap_ex(:,:,hour,1) = Varv_ex;
    end
end

