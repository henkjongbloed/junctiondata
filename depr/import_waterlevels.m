function H = import_waterlevels(RF, DS)

F2 = 'WP2\TwoJunctions\';

addpath(strcat(RF,'Tools\adcptoolsGit'));
addpath(strcat(RF,'Tools\T_Tide'))
addpath(strcat(RF,F2,'\data\processedData\Waterlevels\'));

for i = 1:length(DS)
    H{i,1} = load('NAPWaterLevels1415.mat').(DS{i});
    H{i,1}.name = DS{i};
end
