function [adcp, ctd, h, BN] = import_data(RF, DS, varargin)

% SET PATHS: [SET YOUR OWN]

addpath(strcat(RF,'Tools\adcptoolsGit'));
addpath(strcat(RF,'Tools\loess-master'));
addpath(strcat(RF,'Tools\T_Tide'))
%addpath(strcat(RF,'Tools\example_adcp_processing\muaramuntai25Aug'));
addpath(strcat(RF,'Tools\BrewerMap-master'));
addpath(strcat(RF,'Tools\subaxis'));

F2 = 'WP2\TwoJunctions\';
addpath(strcat(RF,F2,'\data\processedData'));

%addpath(strcat(RF,'processedData'));

% addpath(strcat(RF,'Data\',DS));
% addpath(strcat(RF,'Data\MerelFiles\_processed'));

%addpath(strcat(RF,'DataAnalysis\ProcessedData'));

%addpath(strcat(RF,'Tools\BrewerMap-master'));

%cd(strcat(RF,'DataAnalysis'))

adcp = load(['ADCP',DS], '-mat').adcp;
h = load(['H',DS,'.mat']);
h.time            = datetime(h.date, 'ConvertFrom', 'datenum', 'Format', 'HH:mm:ss');
if strcmp(DS,'NMOMNW15')
    ctd = load(['CTD',DS,'.mat']).CTD_20150914 ;
    BN = {' NM', ' OM', ' NWW'};
else
    ctd = load('OBSsal.mat').OBSsal;
    BN = {' HK', ' OMS', ' OMN'};
end



end
