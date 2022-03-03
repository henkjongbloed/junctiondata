function ctd = hctd(h,ctd,DS)
%Very shady function - modify this.

if strcmp(DS,'NMOMNW15')
    hour              = datenum('01-Jan-2015 01:00:00')-datenum('01-Jan-2015 00:00:00'); %hour = 1/24
    ctd.time = datetime(ctd.Time2+hour, 'ConvertFrom','datenum','Format','HH:mm:ss');
    ctd.time_unique = unique(ctd.time);
    ctd.h          = interp1(h.time, h.level, ctd.time, 'pchip'); %Piecewise Cubic Hermite Interpolating Polynomial. - hoeft niet
    ctd.Z_nap      = - ctd.Z + ctd.h; % correct for the ship height above seabed. DIT CONTROLEREN!!
elseif strcmp(DS, 'OMHA14')
    hour              = datenum('01-Jan-2015 01:00:00')-datenum('01-Jan-2015 00:00:00'); %hour*24 = 1
    ctd.time = datetime(ctd.T+hour, 'ConvertFrom','datenum','Format','HH:mm:ss'); % This is very uncertain!!!
    ctd.time_unique = unique(ctd.time);
    %ctd.Z_nap      = - ctd.Z + ctd.h; % correct for the ship height above seabed. DIT CONTROLEREN!!
    %print('No CTD data found')
    ctd.Trans(1:929) = 1;
    ctd.Trans(930:2639) = 2;
    ctd.Trans(2640:length(ctd.T)) = 3;
    ctd.Hour = getHour(ctd.T, ctd.Trans);
end