function [adcp, ctd, h, dataset_name, tak_names] = import_data2(dataset_idx)
RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'WP2\TwoJunctions\code'])

dataset_name = {'OMHA14','NMOMNW15'}; 
dataset_name = dataset_name{dataset_idx};

tak_names = {{' Hartel Canal'; ' Old Meuse South'; ' Old Meuse North'},...
    {' New Meuse'; ' Old Meuse'; ' Rotterdam Waterway'}}; 

tak_names = tak_names{dataset_idx};
addpath("helpers\", "saved\")
[adcp, ctd, h] = import_data(RF, dataset_name);
end