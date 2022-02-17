function CTD = processSal(V, U, ctd, DS, b)

for hour = 1:13 %Change this
    ind       = find(and(ctd.Hour(:) == hour, ctd.Trans(:) == b)); %Aangepast
    
    %verder aanpassen: Welke CTD data horen bij welke ADCP data??
    
    CTD.X{hour}          =  ctd.X(ind); %LatLon
    CTD.Y{hour}          =  ctd.Y(ind); %LatLon
    CTD.Z{hour}          =  - ctd.Z(ind);% + U.mesh(hour).water_level; %corrected for ship elevation.
    CTD.sal{hour}          =  ctd.S(ind); %salinity
    % WGS2UTM ook bij ADCPData gebruikt?? -> 
    if strcmp(DS,'NMOMNW15')    % Does not work correctly yet
        [CTD.x{hour}, CTD.y{hour}]      = wgs2utm(CTD.Y{hour},CTD.X{hour}); %Convert LatLon to UTM - pay close attention to coordinate systems!!
    elseif strcmp(DS,'OMHA14') % Works correctly. CTD data are in RD coordinates
        [templon, templat, ~]      = rd2wgs(CTD.X{hour},CTD.Y{hour}); %Convert RD to LatLon
        [CTD.x{hour}, CTD.y{hour}]      = wgs2utm(templat, templon); %Convert LatLon to UTM
    end
    
    CTD.P_tmp{hour}      = [CTD.x{hour},CTD.y{hour},CTD.Z{hour}];
    
    % create CTD position data vector
%     cur_mesh = U.mesh(hour);
%     
%     % Transform the CTD data (UTM) to mesh coordinates. Firstly, transform
%     % to the same coordinates as the mesh cells - what are they?
%     %P          =
%     % [ss,ns] = xs.xy2sn(cur_mesh.x_middle(cur_mesh.col_to_cell),cur_mesh.y_middle(cur_mesh.col_to_cell),orig_vel{crp}(:,1),orig_vel{crp}(:,2));
%     [CTD.sctd{hour}, CTD.nctd{hour}] = U.xs.xy2sn( CTD.x{hour}, CTD.y{hour});
%     CTD.zctd{hour} = CTD.Z{hour};
%     
%     CTD.mn{hour} = mean(CTD.nctd{hour});
%     CTD.ms{hour} = mean(CTD.sctd{hour});
%     CTD.velh{hour} = U.vel{hour};
%     CTD.n_cells{hour} = mean(cur_mesh.n_patch,1)';
%     CTD.z_cells{hour} = mean(cur_mesh.z_patch,1)';
%     
%     CTD.Fv{hour}  = scatteredInterpolant(CTD.sctd{hour} - CTD.ms{hour}, CTD.nctd{hour} - CTD.mn{hour}, CTD.zctd{hour}, CTD.sal{hour}, 'linear', 'nearest');
%     CTD.S{hour} = CTD.Fv{hour}(0.*ones(size(CTD.n_cells{hour})),  CTD.n_cells{hour},  CTD.z_cells{hour}); %here, take other coordinates
    
    %eta(hour) = mesh(hour).water_level;
end

end