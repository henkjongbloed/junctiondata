function [h,msh]=importHMSH(y)

if y == '2014'
    % August 2014: Hartelkanaal (Verbeek)
    load H_201408.mat %Imports waterlevels Spijkenisse: 1 Aug 2014 until 31 Aug 2014, 10 minute intervals
    %load ADCP_20140812.mat %1: HAK 2: OMS 3: OMN  , not needed (already imported in MSH_20140812)
    load MSH_20140812.mat %1: HAK 2: OMS 3: OMN (only velocity data)
    h = H_201408;
    msh = MSH_20140812;
elseif y == '2015'
    % September 2015: Nieuwe Waterweg (Poelman)
    load H_2015091415.mat % Imports waterlevels Geulhaven, 14-Sep-15 00:00 until 15-Sep-15 12:24, 1 min interval
    %load CTD_20150914.mat  %See msh_CTD , not needed (already imported in MSH_20150914)
    load MSH_20150914.mat %1: NM 2: OM 3: RWW (salt and velocity data)
    h = H_2015091415;
    msh = MSH_20150914;
elseif y == '2014Test'
     % August 2014: Hartelkanaal (Verbeek)
    load H_201408.mat %Imports waterlevels Spijkenisse: 1 Aug 2014 until 31 Aug 2014, 10 minute intervals
    %load ADCP_20140812.mat %1: HAK 2: OMS 3: OMN  , not needed (already imported in MSH_20140812)
    load MSH_20140812Test.mat %1: HAK 2: OMS 3: OMN (only velocity data)
    h = H_201408;
    msh = MSH_20140812;
elseif y == '2015Test'
       % September 2015: Nieuwe Waterweg (Poelman)
    load H_2015091415.mat % Imports waterlevels Geulhaven, 14-Sep-15 00:00 until 15-Sep-15 12:24, 1 min interval
    %load CTD_20150914.mat  %See msh_CTD , not needed (already imported in MSH_20150914)
    load MSH_20150914Test.mat %1: NM 2: OM 3: RWW (salt and velocity data)
    h = H_2015091415;
    msh = MSH_20150914;
else
    error('Enter correct year')
end
