%Create CTD - THIS SCRIPT IS NOT COMPLETE!!

% load WURData_processed.mat
% load MSH_ADCP_CTD.mat
%% Load the salinity data files
% Identify data files using the common part of the filenames
opt = 'load';
switch opt %Import if seperatate CTD files have to be imported and combined; load if combined file is present.
    case 'import'
        path{1}            = 'C:\Users\jongb013\OneDrive - WageningenUR\PhDHenkJongbloed\2. Programming\DataTools\DataHartelkanaal\14sep2015zout\CTD_NM_proc\'; %Location of CTD data
        path{3}            = 'C:\Users\jongb013\OneDrive - WageningenUR\PhDHenkJongbloed\2. Programming\DataTools\DataHartelkanaal\14sep2015zout\CTD_NW_proc\'; %Location of CTD data
        path{2}            = 'C:\Users\jongb013\OneDrive - WageningenUR\PhDHenkJongbloed\2. Programming\DataTools\DataHartelkanaal\14sep2015zout\CTD_OM_proc\'; %Location of CTD data
        DepName            = 'CC1232014_20150914';            % Common part of the files names
        
        % Create structure to store data in
        CTD.Trans=[];
        CTD.Date=[];
        CTD.Time=[];
        CTD.Time2=[];
        CTD.C=[];
        CTD.T=[];
        CTD.S=[];
        CTD.X=[];
        CTD.Y=[];
        CTD.Z=[];
        CTD.Pressure=[];
        CTD.SC=[];
        CTD.SoundV=[];
        CTD.Density=[];
        CTD.ID=[];
        
        
        for ct=1:3
            addpath(path{ct})
            allFiles=dir([path{ct},DepName,'*']);                  % read all filenames
            allFiles=({allFiles(~[allFiles(:).isdir]).name})'; % We now have a vector containing all filenames.
            
            %    Loop to load all files and store data in cell arrays
            for i=1:size(allFiles)
                load(allFiles{i});
                CTD.Trans=[CTD.Trans; repmat(ct, size(Temperature))];
                CTD.Date=[CTD.Date; repmat(str2num(datestr(CastTimeUtc', 'YYYYMMDD')), size(Temperature))];
                CTD.Time2=[CTD.Time2; repmat(datenum(CastTimeUtc'), size(Temperature))];
                CTD.Time=[CTD.Time; repmat(str2num(datestr(CastTimeUtc', 'HHMMSS')), size(Temperature))];
                CTD.C=[CTD.C; Conductivity];
                CTD.T=[CTD.T; Temperature];
                CTD.S=[CTD.S; Salinity];
                CTD.X=[CTD.X; repmat(nanmean([LongitudeEnd, LongitudeStart]), size(Temperature))];
                CTD.Y=[CTD.Y; repmat(nanmean([LatitudeEnd, LatitudeStart]), size(Temperature))];
                CTD.Z=[CTD.Z; Depth];
                CTD.Pressure=[CTD.Pressure; Pressure];
                CTD.SC=[CTD.SC; Specific_conductance];
                CTD.SoundV=[CTD.SoundV; Sound_velocity];
                CTD.Density=[CTD.Density; Density];
                CTD.ID=[CTD.ID; repmat(i, size(Temperature))];
                %     cellsz = structfun(@size,CTD,'uni',false);
                %     Size=[Size,
                %     clearvars -except adcp allFiles CTD DepName path tran i msh UTM ct%Clear all imported data
            end
            rmpath(path{ct});    %Remove path when all files stored in a folder are imported, so there is no risk of overwriting.
        end
        
        % Assign values for hours (1:13) based on time.
        hours = [50800; 61100; 70300; 75000; 85500; 100000; 105400; 115000; 130000; 140000; 145100; 161000; 170000; 180000]; %look into this
        CTD.Hour=zeros(length(CTD.Time),1);
        k = 1:13;
        for i = 1:length(CTD.Time)
            dat = CTD.Time(i);
            if size(CTD.Time,1)>1
                dat = CTD.Time(i); % dataT(1,i)
            end
            el  = find( hours(k+1) > dat & hours(k) <= dat );
            if el
                CTD.Hour(i) = el; %Should this be k?
            end
            
        end
        
        
        % [CTD.UTM(:,1),CTD.UTM(:,2)] = wgs2utm(CTD.Y,CTD.X);
        % [CTD.UTM(:,1),CTD.UTM(:,2)]=ll2utm(CTD.Y,CTD.X, 'int24');
        save('CTD_v3', 'CTD')
        %
    case 'load'
        load CTD_v3.mat
end
