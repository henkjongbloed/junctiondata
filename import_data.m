function [adcp, ctd, h] = import_data(RF, DS)

addpath(strcat(RF,'Tools\adcptoolsGit'));
addpath(strcat(RF,'Tools\loess-master'));
addpath(strcat(RF,'Tools\T_Tide'))
addpath(strcat(RF,'Tools\BrewerMap-master'));
addpath(strcat(RF,'Tools\subaxis'));

F2 = 'WP2\TwoJunctions\';
addpath(strcat(RF,F2,'\data\processedData'));

adcp = load(['ADCP',DS], '-mat').adcp;
h = load(['H',DS,'.mat']);
ctd = load(['CTD', DS], '-mat').ctd;
h.time            = datetime(h.date, 'ConvertFrom', 'datenum', 'Format', 'HH:mm:ss');
% if strcmp(DS,'NMOMNW15')
%     ctd = load(['CTD',DS,'.mat']).CTD_20150914 ;
%     BN = {' NM'; ' OM'; ' NWW'};
% else
%     ctd = load('OBSsal.mat').OBSsal;
%     BN = {' HK'; ' OMS'; ' OMN'};
% end
% 
% ctd = split_branches(ctd, DS, BN);
% ctd = split_branches(ctd, DS, BN);
% save(string(strcat('CTD',DS)) , 'ctd')

end
